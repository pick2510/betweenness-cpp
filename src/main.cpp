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
#include "natural_sort.hpp"
#include "graph.h"
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/status.hpp>

constexpr int MASTER = 0;

int main(int argc, char **argv)
{
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  std::vector<Result> results;
  Config runningConf;
  std::map<int, std::string> inv_vertice_map;

  if (world.rank() == MASTER)
  {
    
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
    auto chain_file_list = glob_deq(chainpattern);
    std::vector<std::string> radius_file_list = glob(xyzpattern);

    // sanity check of file list

    if (radius_file_list.empty())
    {
      BOOST_LOG_TRIVIAL(info) << "Need at least one file with radii *.tet\n";
      exit(EXIT_FAILURE);
    }
    SI::natural::sort(chain_file_list);

    // Read one radius file, create lookup table for continous particle ids.
    // Create
    auto radius_file = read_file(radius_file_list[0]);
    radius_file.pop_back();
    std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
    auto lookup_table = get_lookup_table(radius_file);
    auto vertice_map = get_vertice_map(radius_file);
    auto inv_vertice_map = inverse_map(vertice_map);
    auto t_len = chain_file_list.size();
    long v_index = 0;
   
    //std::atomic<long> index{0};
    //double percent{0.0};
    results.resize(t_len);
    std::vector<boost::mpi::request> reqs(world.size());

    for (int dst_rank = 1; dst_rank < world.size(); ++dst_rank)
    {
      Job job{chain_file_list.front(), vertice_map};
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job "
                              << " to SLAVE " << dst_rank << "\n";
      world.isend(dst_rank, 0, job);
      chain_file_list.pop_front();
      // Post receive request for new jobs requests by slave [nonblocking]
      reqs[dst_rank] = world.irecv(dst_rank, 0, results[v_index++]);
    }

    while (chain_file_list.size() > 0)
    {
      bool stop;
      for (int dst_rank = 1; dst_rank < world.size(); ++dst_rank)
      {
        // Check if dst_rank is done
        if (reqs[dst_rank].test())
        {
          BOOST_LOG_TRIVIAL(info) << "[MASTER] Rank " << dst_rank << " is done.\n";
          // Check if there is remaining jobs
          if (chain_file_list.size() > 0)
          {
            // Tell the slave that a new job is coming.
            Job job{chain_file_list.front(), vertice_map};
            stop = false;
            world.isend(dst_rank, 0, stop);
            // Send the new job.
            BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending new job ("
                      << ") to SLAVE " << dst_rank << ".\n";
            world.isend(dst_rank, 0, job);
            chain_file_list.pop_front();
            reqs[dst_rank] = world.irecv(dst_rank, 0, results[v_index++]);
          }
          else
          {
            // Send stop message to slave.
            stop = true;
            world.isend(dst_rank, 0, stop);
          }
        }
      }
      usleep(1000);
    }
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Sent all jobs.\n";

     // Listen for the remaining jobs, and send stop messages on completion.
    bool all_done = false;
    while (!all_done) {
      all_done = true;
      for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
        if (reqs[dst_rank].test()) {
            // Tell the slave that it can exit.
            bool stop = true;
            world.isend(dst_rank, 0, stop);
        }
        else {
          all_done = false;
        }
      }
      usleep(1000);
    }
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Handled all jobs, killed every process.\n";


    /*
  for (std::size_t i = 0; i < chain_file_list.size(); ++i)
  {
    Graph mygraph(chain_file_list[i], vertice_map);
    mygraph.calc();

    results.push_back(mygraph.get_result());
    index++;
    percent = (index / t_len) * 100;
    BOOST_LOG_TRIVIAL(info) << percent << "% (" << index << " of " << t_len << ") done";
  }
  */
  } else {

    bool stop = false;
    while(!stop) {
      // Wait for new job
      Job job;
      world.recv(0, 0, job);
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << world.rank()
                << "] Received job " << " from MASTER.\n";
      // Perform "job"
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: "<< world.rank() 
                << "Start Job\n";
      Graph mygraph(job.file, job.v_map);
      mygraph.calc();
      auto res = mygraph.get_result();
      // Notify master that the job is done
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << world.rank() 
                << "] Done with job " << ". Notifying MASTER.\n";
      world.send(0, 0, res);
      // Check if a new job is coming
      world.recv(0, 0, stop);
    }

  BOOST_LOG_TRIVIAL(info) << "~~~~~~~~ Rank " << world.rank() << " is exiting ~~~~~~~~~~~\n";


  }
  if (world.rank() == MASTER){
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
