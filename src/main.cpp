#include <iostream>
#include <iomanip>
#include <string>
#include <deque>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <atomic>

#include "main.h"
#include "utils.h"
#include "data.h"
#include "natural_sort.hpp"
#include "graph.h"
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/status.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

constexpr int MASTER = 0;
constexpr int TAG_RESULT = 1;
constexpr int TAG_BREAK = 2;
constexpr int TAG_FILE = 10;


int main(int argc, char **argv)
{
  //Initialize global object
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  std::vector<Result> results;
  Config runningConf;
  std::map<int, std::string> inv_vertice_map;
  std::map<std::string, int> vertice_map;
  std::vector<std::string> keys;
  std::vector<int> vals;
  std::vector<std::string> radius_file_list;
  std::deque<std::string> chain_file_list;

  if (world.rank() == MASTER)
  {
    //( )
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Node betweenness centrality";
    BOOST_LOG_TRIVIAL(info) << "using BOOST Graph Library";
    BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
    BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    try
    {
      runningConf = getCL(argc, argv);
    }
    catch (std::invalid_argument e)
    {
      BOOST_LOG_TRIVIAL(info) << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }

    // define pattern fro globbing chain files, trim strings, and glob

    std::string chainpattern(runningConf.InputPath + "/*.chain");
    std::string xyzpattern(runningConf.InputPath + "/*.tet");
    xyzpattern = trim(xyzpattern);
    chainpattern = trim(chainpattern);
    chain_file_list = glob_deq(chainpattern);
    radius_file_list = glob(xyzpattern);
    auto radius_file = read_file(radius_file_list[0]);
    radius_file.pop_back();
    std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
    auto lookup_table = get_lookup_table(radius_file);
    vertice_map = get_vertice_map(radius_file);
    inv_vertice_map = inverse_map(vertice_map);
    keys = getKeys(vertice_map);
    vals = getVals(vertice_map);
  }

  //Broadcast keys and vals vector for vertice_map reconstruction.
  broadcast(world, keys, MASTER);
  broadcast(world, vals, MASTER);

  if (world.rank() == MASTER)
  {
    // sanity check of file list

    if (radius_file_list.empty())
    {
      BOOST_LOG_TRIVIAL(info) << "Need at least one file with radii *.tet\n";
      exit(EXIT_FAILURE);
    }
    SI::natural::sort(chain_file_list);

    auto t_len = chain_file_list.size();
    //Check that processes are <= length of tasks
    if (world.size() > t_len + 1)
    {
      BOOST_LOG_TRIVIAL(error) << "Too many processes (" << world.size()
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use " << t_len + 1 << " ranks or less\n";
      return 0;
    }
    long v_index = 0;
    results.resize(t_len);
    std::vector<boost::mpi::request> reqs(world.size());

    for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank)
    {
      std::string file{chain_file_list.front()};
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job "
                              << file << " to SLAVE (first loop) " << dst_rank << "\n";

      world.isend(dst_rank, TAG_FILE, file);
      chain_file_list.pop_front();
      // Post receive request for new jobs requests by slave [nonblocking]
      reqs[dst_rank] = world.irecv(dst_rank, TAG_RESULT, results[v_index++]);
    }

    while (chain_file_list.size() > 0)
    {
      bool stop;
      for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank)
      {
        // Check if dst_rank is done
        if (reqs[dst_rank].test())
        {
          BOOST_LOG_TRIVIAL(info) << "[MASTER] Rank " << dst_rank << " is done.\n";
          // Check if there is remaining jobs
          if (chain_file_list.size() > 0)
          {
            // Tell the slave that a new job is coming.
            stop = false;
            world.isend(dst_rank, TAG_BREAK, stop);
            std::string file{chain_file_list.front()};
            // Send the new job.
            BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending new job ("
                                    << file << ") to SLAVE " << dst_rank << ".\n";

            world.isend(dst_rank, TAG_FILE, file);
            chain_file_list.pop_front();
            reqs[dst_rank] = world.irecv(dst_rank, TAG_RESULT, results[v_index++]);
          }
          else
          {
            // Send stop message to slave.
            stop = true;
            world.isend(dst_rank, TAG_BREAK, stop);
          }
        }
      }
      usleep(1000);
    }
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Sent all jobs.\n";

    // Listen for the remaining jobs, and send stop messages on completion.
    bool all_done = false;
    while (!all_done)
    {
      all_done = true;
      for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank)
      {
        if (reqs[dst_rank].test())
        {
          // Tell the slave that it can exit.
          bool stop = true;
          world.isend(dst_rank, TAG_BREAK, stop);
        }
        else
        {
          all_done = false;
        }
      }
      usleep(1000);
    }
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Handled all jobs, killed every process.\n";
  }
  // SLAVE CODE
  if (world.rank() != MASTER)
  {
    //Generate Map from string and int vector
    std::map<std::string, int> v_map;
    std::transform(keys.begin(), keys.end(), vals.begin(), std::inserter(v_map, v_map.end()), [](std::string a, int b) {
      return std::make_pair(a, b);
    });
    bool stop = false;
    while (!stop)
    {
      std::string file;

      BOOST_LOG_TRIVIAL(info) << "Initialized Job";
      world.recv(0, TAG_FILE, file);

      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << world.rank()
                              << "] Received job "
                              << file << " from MASTER.\n";
      // Perform "job"
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << world.rank()
                              << "] Start Job\n";
      Graph mygraph(file, v_map);
      mygraph.calc();
      auto res = mygraph.get_result();
      // Send result
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << world.rank()
                              << "] Done with job "
                              << file << ". Send Result.\n";
      world.send(0, TAG_RESULT, res);
      // Check if a new job is coming
      world.recv(0, TAG_BREAK, stop);
    }

    BOOST_LOG_TRIVIAL(info) << "Rank " << world.rank() << " is exiting\n";
  }
  if (world.rank() == MASTER)
  {
    std::sort(results.begin(), results.end(), cmp_ts);
    std::ofstream ts_mean_file(runningConf.OutputPath + "/properties.csv");
    write_ts_header(ts_mean_file, runningConf);
    ts_mean_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (auto &v : results)
    {
      ts_mean_file << std::to_string(v.ts) << runningConf.sep << std::to_string(v.mean)
                   << runningConf.sep << std::to_string(v.var) << runningConf.sep << std::to_string(std::sqrt(v.var))
                   << runningConf.sep << std::to_string(v.skew) << runningConf.sep << std::to_string(v.kur)
                   << runningConf.sep << std::to_string(v.q_090) << runningConf.sep << std::to_string(v.q_099)
                   << "\n";
      std::ofstream ts_file(runningConf.OutputPath + "/centrality_" + std::to_string(v.ts) + ".csv");
      write_cent_header(ts_file, runningConf);
      ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
      for (auto &kv : v.b_centrality)
      {
        ts_file << inv_vertice_map[kv.first] << runningConf.sep << kv.second << "\n";
      }
      ts_file.close();
    }
    ts_mean_file.flush();
    ts_mean_file.close();
  }
}
