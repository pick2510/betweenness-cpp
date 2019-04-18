#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "data.h"
#include "graph.h"
#include "main.h"
#include "natural_sort.hpp"
#include "utils.h"
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/status.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <Eigen/Eigen>
#include <unistd.h>

constexpr int MASTER = 0;
constexpr int TAG_RESULT = 1;
constexpr int TAG_BREAK = 2;
constexpr int TAG_FILE = 10;
constexpr int HOSTNAME_LEN = 255;
constexpr int TAG_PART_TS = 20;

int main(int argc, char **argv)
{
  // Initialize global object
  char hostname[HOSTNAME_LEN]{};
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  Config runningConf{};
  std::map<int, std::string> inv_vertice_map{};
  std::map<std::string, int> vertice_map{};
  std::vector<std::string> keys{};
  std::vector<int> vals {};
  std::vector<Result> results {};
  std::vector<std::string> radius_file_list {};
  std::deque<std::string> chain_file_list {};
  double *ts_particle;
  auto rank = world.rank();
  auto world_size = world.size();
  int t_len{};
  gethostname(hostname, HOSTNAME_LEN);

  if (rank == MASTER) {
    // MASTER CODE
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Node betweenness centrality";
    BOOST_LOG_TRIVIAL(info) << "using BOOST Graph Library";
    BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
    BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
    try {
      runningConf = getCL(argc, argv);
    }
    catch (std::invalid_argument e) {
      BOOST_LOG_TRIVIAL(info) << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }

    // define pattern fro globbing chain files, trim strings, and glob

    std::string chainpattern(runningConf.InputPath + "/postchain/*.chain");
    std::string xyzpattern(runningConf.InputPath + "/postxyz/*.tet");
    xyzpattern = trim(xyzpattern);
    chainpattern = trim(chainpattern);
    chain_file_list = glob_deq(chainpattern);
    radius_file_list = glob(xyzpattern);
    if (radius_file_list.empty()) {
      BOOST_LOG_TRIVIAL(info) << "Need at least one file with radii *.tet\n";
      world.abort(-255);
      exit(EXIT_FAILURE);
    }
    auto radius_file = read_file(radius_file_list[0]);

    radius_file.pop_back();
    std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
    auto lookup_table = get_lookup_table(radius_file);
    vertice_map = get_vertice_map(radius_file);
    inv_vertice_map = inverse_map(vertice_map);
    keys = getKeys(vertice_map);
    vals = getVals(vertice_map);
  }

  // Broadcast keys and vals vector for vertice_map reconstruction.
  broadcast(world, keys, MASTER);
  broadcast(world, vals, MASTER);


  //MASTER CODE

  if (rank == MASTER) {
    // sanity check of file list
    SI::natural::sort(chain_file_list);
    t_len = chain_file_list.size();
    BOOST_LOG_TRIVIAL(info) << "t_len = " << t_len;
    results.resize(chain_file_list.size());
    auto p_size = vertice_map.size();
    // Check that processes are <= length of tasks
    if (world_size > t_len + 1) {
      BOOST_LOG_TRIVIAL(error) << "Too many processes (" << world_size
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use " << t_len + 1 << " ranks or less\n";
      return 0;
    }
    ts_particle = new double[t_len * p_size] {};
    BOOST_LOG_TRIVIAL(info)
        << "Size of ts_particle: " << p_size * t_len * sizeof(double);
    long v_index = 0;
    BOOST_LOG_TRIVIAL(info) << "Memory: " << results.size();
    std::vector<boost::mpi::request> reqs_world(world_size);
    std::vector<boost::mpi::request> reqs_ts(reqs_world);

    for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
      std::string file{chain_file_list.front()};
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job " << file
                              << " to SLAVE (first loop) " << dst_rank << "\n";
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (v_index / t_len) * 100.0 << "% done";
      world.send(dst_rank, TAG_FILE, file.data(), file.size());
      chain_file_list.pop_front();
      // Post receive request for new jobs requests by slave [nonblocking]
      reqs_world[dst_rank] =
          world.irecv(dst_rank, TAG_RESULT, results[v_index]);
      reqs_ts[dst_rank] = world.irecv(dst_rank, TAG_PART_TS,
                                      &ts_particle[v_index++ * p_size], p_size);
    }
    bool stop = false;
    while (chain_file_list.size() > 0) {
      auto status = wait_any(reqs_world.begin(), reqs_world.end());
      auto dst_rank = status.first.source();
      reqs_ts[dst_rank].wait();
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Rank " << dst_rank << " is done.\n";
      // Check if there is remaining jobs

      // Tell the slave that a new job is coming.
      world.send(dst_rank, TAG_BREAK, stop);

      // Send the new job.
      std::string file{chain_file_list.front()};
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] Sending new job (" << file << ") to SLAVE " << dst_rank;
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (v_index / t_len) * 100.0 << "% done";
      world.send(dst_rank, TAG_FILE, file.data(), file.size());
      chain_file_list.pop_front();
      reqs_world[dst_rank] =
          world.irecv(dst_rank, TAG_RESULT, results[v_index]);
      reqs_ts[dst_rank] = world.irecv(dst_rank, TAG_PART_TS,
                                      &ts_particle[v_index++ * p_size], p_size);
    }
   
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Sent all jobs.\n";
    wait_all(reqs_world.begin(), reqs_world.end());
    wait_all(reqs_ts.begin(), reqs_ts.end());
    stop = true;
    for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
      world.send(dst_rank, TAG_BREAK, stop);
    }

    BOOST_LOG_TRIVIAL(info)
        << "[MASTER] Handled all jobs, killed every process.\n";
  }
  
  
  
  // SLAVE CODE
  if (rank > MASTER) {
    // Generate Map from string and int vector
    std::map<std::string, int> v_map = constructMap(keys, vals);
    bool stop = false;
    while (!stop) {
      auto status = world.probe(0, TAG_FILE);
      auto nbytes = status.count<char>();
      char *f_recv = new char[nbytes.get() + 1]{};
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Intialized JOB";
      world.recv(0, TAG_FILE, f_recv, nbytes.get());
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << rank << "] (" << hostname
                              << ") Recevied: " << nbytes.get() << " bytes";
      std::string file(f_recv);
      delete[] f_recv;
      BOOST_LOG_TRIVIAL(info) << file;

      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << rank << "] (" << hostname
                              << ") Received job " << file << " from MASTER.\n";
      // Perform "job"
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Start Job\n";
      Graph mygraph(file, v_map);
      mygraph.calc();
      auto res = mygraph.get_result();
      // Send result
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Done with job "
          << file << ". Send Result.\n";
      world.send(0, TAG_RESULT, res);
      world.send(0, TAG_PART_TS, res.vals.data(), res.vals.size());

      // Check if a new job is coming
      world.recv(0, TAG_BREAK, stop);
    }

    BOOST_LOG_TRIVIAL(info) << "Rank " << rank << " is exiting\n";
  }
  
  //MASTTER CODE
  
  if (rank == MASTER) {
    std::sort(results.begin(), results.end(), cmp_ts);
    std::ofstream ts_mean_file(runningConf.OutputPath + "/properties.csv");
    write_ts_header(ts_mean_file, runningConf);
    ts_mean_file << std::setprecision(std::numeric_limits<double>::digits10 +
                                      1);
    output_ts(ts_mean_file, runningConf, results, inv_vertice_map);
    ts_mean_file.flush();
    ts_mean_file.close();
    BOOST_LOG_TRIVIAL(info) << "TEST PARTICLE:  " << ts_particle[12];
    Eigen::Map<Eigen::MatrixXd> particle_matrix(ts_particle, t_len, vertice_map.size());

    delete[] ts_particle;
   
  }
}
