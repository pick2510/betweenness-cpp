

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

#include "INIReader.h"
#include "betweenness_sqlite.h"
#include "data.h"
#include "graph_sqlite.h"
#include "grid.h"
#include "memory"
#include "natural_sort.hpp"
#include "utils.h"
#include <Eigen/Eigen>
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
#include <unistd.h>

int main(int argc, char **argv)
{
  // Initialize global object
  char hostname[HOSTNAME_LEN]{};
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  boost::mpi::communicator tscom(world, boost::mpi::comm_duplicate);
  Config runningConf{};
  std::string configPath;
  std::map<int, double> radius_map;
  std::map<int, std::string> inv_vertice_map{};
  std::map<std::string, int> vertice_map{};
  std::vector<std::string> keys{};
  std::vector<int> vals{};
  std::vector<Betweenness_Result> results{};
  std::deque<long> ts{};
  std::vector<long> ts_res{};
  std::string path_contact{};
  std::string path_particle{};
  std::vector<decomp_table> decomp_tab{};
  std::unique_ptr<double[]> ts_particle;
  auto rank = world.rank();
  auto world_size = world.size();
  int t_len{};
  gethostname(hostname, HOSTNAME_LEN);

  if (rank == MASTER) {
    INIReader reader;
    // MASTER CODE
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Node betweenness centrality";
    BOOST_LOG_TRIVIAL(info) << "using BOOST Graph Library";
    BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
    BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
    try {
      configPath = getConfigPath(argc, argv);
    }
    catch (std::invalid_argument &e) {
      BOOST_LOG_TRIVIAL(info) << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    try {
      reader = parseConfigFile(configPath);
    }
    catch (std::invalid_argument &e) {
      BOOST_LOG_TRIVIAL(info) << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    runningConf = getBetweennessConfigObj(reader, "betweenness_sqlite");
    // define pattern fro globbing chain files, trim strings, and glob

    path_contact =
        std::string{runningConf.InputPath + "/" + runningConf.contact_filename};
    path_particle = std::string{runningConf.InputPath + "/" +
                                runningConf.particle_filename};
    auto ts_store = initTSStorage(path_contact);
    auto contact_store = indexContactStorage(path_contact);
    auto particle_store = ParticleIndexStorage(path_particle);
    auto decomp_db = initDecompstorage(path_contact);
    decomp_tab = decomp_db.get_all<decomp_table>();
    contact_store.sync_schema();
    particle_store.sync_schema();
    auto radius_storage = initRadstorage(path_contact);
    for (auto &elem : radius_storage.iterate<radius>()) {
      radius_map[elem.id] = elem.rad;
    }

    auto ts_list = ts_store.get_all<ts_column>();
    for (const auto &elem : ts_list) {
      ts.push_back(elem.ts);
    }
    std::sort(ts.begin(), ts.end());
    auto lookup_table = get_lookup_table(radius_map);
    vertice_map = get_vertice_map(radius_map);
    inv_vertice_map = inverse_map(vertice_map);
    keys = getKeys(vertice_map);
    vals = getVals(vertice_map);
  }

  // Broadcast keys and vals vector for vertice_map reconstruction.
  broadcast(world, keys, MASTER);
  broadcast(world, vals, MASTER);
  broadcast(world, decomp_tab, MASTER);

  // MASTER CODE

  if (rank == MASTER) {
    // sanity check of file list
    t_len = ts.size();
    BOOST_LOG_TRIVIAL(info) << "t_len = " << t_len;
    results.resize(t_len);
    auto p_size = vertice_map.size();
    // Check that processes are <= length of tasks
    if (world_size > t_len + 1) {
      BOOST_LOG_TRIVIAL(error) << "Too many processes (" << world_size
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use " << t_len + 1 << " ranks or less\n";
      world.abort(EXIT_FAILURE);
      return 0;
    }
    if (world_size < 2) {
      BOOST_LOG_TRIVIAL(error) << "Too few processes (" << world_size
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use at least 2 ranks\n";

      world.abort(EXIT_FAILURE);
      return 0;
    }
    ts_particle = std::make_unique<double[]>(t_len * p_size);
    BOOST_LOG_TRIVIAL(info)
        << "Size of ts_particle: " << p_size * t_len * sizeof(double);
    long v_index = 0;
    BOOST_LOG_TRIVIAL(info) << "Memory: " << results.size();
    std::vector<boost::mpi::request> reqs_world(world_size);
    std::vector<boost::mpi::request> reqs_ts(reqs_world);

    for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
      long timestep{ts.front()};
      auto contact_db = indexContactStorage(path_contact);
      auto particle_db = ParticleIndexStorage(path_particle);
      auto contact_columns = contact_db.get_all<ContactColumns>(
          where(c(&ContactColumns::ts) == timestep));
      auto particle_columns = particle_db.get_all<ParticleColumns>(
          where(c(&ParticleColumns::ts) == timestep));
      auto contact_col_size = contact_columns.size();
      auto particle_col_size = particle_columns.size();
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job " << timestep
                              << " to SLAVE (first loop) " << dst_rank << "\n";
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (v_index / t_len) * 100.0 << "% done";
      world.send(dst_rank, TAG_TS, timestep);
      world.send(dst_rank, TAG_SIZE, contact_col_size);
      world.send(dst_rank, TAG_FILE, contact_columns);
      world.send(dst_rank, TAG_SIZE, particle_col_size);
      world.send(dst_rank, TAG_P_COL, particle_columns);
      ts.pop_front();
      // Post receive request for new jobs requests by slave [nonblocking]
      reqs_world[dst_rank - 1] =
          world.irecv(dst_rank, TAG_RESULT, results[v_index]);
      reqs_ts[dst_rank - 1] = tscom.irecv(
          dst_rank, TAG_PART_TS, &ts_particle[v_index++ * p_size], p_size);
    }
    bool stop = false;
    while (ts.size() > 0) {
      for (unsigned int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
        if (reqs_world[dst_rank - 1].test()) {
          reqs_ts[dst_rank - 1].wait();
          BOOST_LOG_TRIVIAL(info)
              << "[MASTER] Rank " << dst_rank << " is done.\n";
          // Check if there is remaining jobs

          // Tell the slave that a new job is coming.
          world.send(dst_rank, TAG_BREAK, stop);

          // Send the new job.
          long timestep{ts.front()};
          auto contact_db = indexContactStorage(path_contact);
          auto particle_db = ParticleIndexStorage(path_particle);
          auto contact_columns = contact_db.get_all<ContactColumns>(
              where(c(&ContactColumns::ts) == timestep));
          auto particle_columns = particle_db.get_all<ParticleColumns>(
              where(c(&ParticleColumns::ts) == timestep));
          auto contact_col_size = contact_columns.size();
          auto particle_col_size = particle_columns.size();
          BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending new job (" << timestep
                                  << ") to SLAVE " << dst_rank;
          BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
          BOOST_LOG_TRIVIAL(info)
              << "[MASTER] " << (v_index / t_len) * 100.0 << "% done";
          world.send(dst_rank, TAG_TS, timestep);
          world.send(dst_rank, TAG_SIZE, contact_col_size);
          world.send(dst_rank, TAG_FILE, contact_columns);
          world.send(dst_rank, TAG_SIZE, particle_col_size);
          world.send(dst_rank, TAG_PART_TS, particle_col_size);
          ts.pop_front();
          reqs_world[dst_rank - 1] =
              world.irecv(dst_rank, TAG_RESULT, results[v_index]);
          reqs_ts[dst_rank - 1] = tscom.irecv(
              dst_rank, TAG_PART_TS, &ts_particle[v_index++ * p_size], p_size);
        }
      }
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
      long ts;
      std::vector<ContactColumns> contact_cols{};
      std::vector<ParticleColumns> particle_cols{};
      size_t contact_col_size;
      size_t particle_col_size;
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Intialized JOB";
      world.recv(0, TAG_TS, ts);
      BOOST_LOG_TRIVIAL(info) << ts;
      world.recv(0, TAG_SIZE, contact_col_size);
      contact_cols.resize(contact_col_size);
      world.recv(0, TAG_FILE, contact_cols);
      world.recv(0, TAG_SIZE, particle_col_size);
      particle_cols.resize(particle_col_size);
      world.recv(0, TAG_P_COL, particle_cols);

      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << rank << "] (" << hostname
                              << ") Received job " << ts << " from MASTER.\n";
      // Perform "job"
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Start Job\n";
      GraphSQLite mygraph(v_map, ts, contact_cols, decomp_tab, particle_cols);
      mygraph.calc();
      auto res = mygraph.get_result();
      // Send result
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Done with job " << ts
          << ". Send Betweenness_Result.\n";
      world.send(0, TAG_RESULT, res);
      tscom.send(0, TAG_PART_TS, res.vals.data(), res.vals.size());
      // Check if a new job is coming
      world.recv(0, TAG_BREAK, stop);
    }

    BOOST_LOG_TRIVIAL(info) << "Rank " << rank << " is exiting\n";
  }

  // MASTTER CODE

  if (rank == MASTER) {
    std::sort(results.begin(), results.end(), cmp_ts);
    ts_res.reserve(results.size());
    std::for_each(
        results.begin(), results.end(),
        [&ts_res](Betweenness_Result const &res) { ts_res.push_back(res.ts); });
    std::ofstream ts_mean_file(runningConf.OutputPath +
                               "/betweenness_global.csv");
    write_ts_header(ts_mean_file, runningConf);
    ts_mean_file << std::setprecision(std::numeric_limits<double>::digits10 +
                                      1);
    output_centrality_ts(ts_mean_file, runningConf, results, inv_vertice_map);
    ts_mean_file.close();
    BOOST_LOG_TRIVIAL(info) << "TEST PARTICLE:  " << ts_particle[12];
    Eigen::Map<Eigen::MatrixXd> particle_matrix(ts_particle.get(), t_len,
                                                vertice_map.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp single
    {
#pragma omp task
#endif
      output_particle_ts(runningConf, particle_matrix, inv_vertice_map, ts_res);
#if defined(_OPENMP)
#pragma omp task
#endif
      output_cell_ts(runningConf, results, ts_res, decomp_tab);
#if defined(_OPENMP)
    }
#endif
  }
}
