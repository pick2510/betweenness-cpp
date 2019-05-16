//
// Created by std on 10.05.19.
//

#include "properties.h"
#include "INIReader.h"
#include "PotentialEnergy.h"
#include "accumulator.h"
#include "data.h"
#include "grid.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/status.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <map>
#include <serialize_tuple.h>
#include <tuple>

using namespace sqlite_orm;

int main(int argc, char *argv[])
{
  char hostname[HOSTNAME_LEN]{};
  boost::mpi::environment env{argc, argv};
  boost::mpi::communicator world;
  Config runningConf{};
  std::string configPath;
  std::map<int, double> radius_map;
  auto rank = world.rank();
  auto world_size = world.size();
  int t_len{};
  std::deque<long> ts{};
  std::string path;
  decom_vec_storage_t decomp_str;
  std::vector<int> keys{};
  std::vector<double> vals{};
  std::vector<aggr_result_t> results{};
  gethostname(hostname, HOSTNAME_LEN);
  boost::filesystem::path system_path{};
  boost::filesystem::path cellstr_path{};

  if (rank == MASTER) {
    INIReader reader;
    // MASTER CODE
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Potential Energy Calculator";
    BOOST_LOG_TRIVIAL(info) << "for DEM Output";
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
    runningConf = getPropertiesConfigObj(reader);

    path =
        std::string{runningConf.InputPath + "/" + runningConf.contact_filename};
    cellstr_path = boost::filesystem::path(runningConf.OutputPath + "/cells");
    system_path = boost::filesystem::path(runningConf.OutputPath + "/system");
    check_path(cellstr_path);
    check_path(system_path);
    auto stor = initTSStorage(path);
    auto particles = indexContactStorage(path);
    auto radius_storage = initRadstorage(path);
    for (auto &elem : radius_storage.iterate<radius>()) {
      radius_map[elem.id] = elem.rad;
    }
    auto ts_list = stor.get_all<ts_column>();
    for (const auto &elem : ts_list) {
      ts.push_back(elem.ts);
    }

    std::sort(ts.begin(), ts.end());
    auto upper =
        std::lower_bound(ts.begin(), ts.end(), runningConf.spinup_time);
    if (upper != ts.end())
      ts.erase(ts.begin(), upper);
    auto lookup_table = get_lookup_table(radius_map);
    auto cell_db = initDecompstorage(path);
    decomp_str = cell_db.get_all<decomp_table>();
    keys = getKeys(radius_map);
    vals = getVals(radius_map);
  }

  broadcast(world, keys, MASTER);
  broadcast(world, decomp_str, MASTER);
  broadcast(world, vals, MASTER);
  broadcast(world, path, MASTER);
  auto db = indexContactStorage(path);

  if (rank == MASTER) {
    auto db = indexContactStorage(path);
    long counter = 0;
    // sanity check of file list
    t_len = ts.size();
    int chunk_len =
        runningConf.chunk_len > t_len ? t_len : runningConf.chunk_len;
    BOOST_LOG_TRIVIAL(info) << "t_len = " << t_len;
    BOOST_LOG_TRIVIAL(info) << "chunk_len = " << chunk_len;
    results.resize(chunk_len);
    // Check that processes are <= length of tasks
    if (world_size > t_len + 1) {
      BOOST_LOG_TRIVIAL(error) << "Too many processes (" << world_size
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use " << t_len + 1 << " ranks or less\n";
      world.abort(EXIT_FAILURE);
      return EXIT_FAILURE;
    }
    if (world_size < 2) {
      BOOST_LOG_TRIVIAL(error) << "Too few processes (" << world_size
                               << ") for the number of jobs!\n";
      BOOST_LOG_TRIVIAL(error) << "Use at least 2 ranks\n";

      world.abort(EXIT_FAILURE);
      return EXIT_FAILURE;
    }
    long v_index = 0;
    BOOST_LOG_TRIVIAL(info) << "Memory: " << results.size();
    std::vector<boost::mpi::request> reqs_world;
    BOOST_LOG_TRIVIAL(info) << "Intitalize systemfile";
    initialize_output_files(runningConf, decomp_str, system_path, cellstr_path);
    submitCompleteWorld(world, world_size, t_len, ts, db, reqs_world, results,
                        v_index, counter);
    bool stop = false;

    if (!ts.empty()) {
      submitPieces(world, runningConf, t_len, ts, decomp_str, results,
                   system_path, cellstr_path, db, chunk_len, v_index,
                   reqs_world, stop, world_size, counter);
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sent all jobs.\n";
      wait_all(reqs_world.begin(), reqs_world.end());
    }
    else {
      wait_all(reqs_world.begin(), reqs_world.end());
    }
    write_results(runningConf, decomp_str, results, system_path, cellstr_path);
    stop = true;
    for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
      world.send(dst_rank, TAG_BREAK, stop);
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] Handled all jobs, killed every process.\n";
    }
  }
  if (rank > MASTER) {
    std::map<int, double> radius = constructMap(keys, vals);
    bool stop = false;
    while (!stop) {
      tuple_storage_t data_recv{};
      size_t col_size;
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Intialized JOB";
      world.recv(0, TAG_SIZE, col_size);
      data_recv.resize(col_size);
      world.recv(0, TAG_FILE, data_recv);
      auto timestep = std::get<10>(*data_recv.begin());
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] recevied ts: " << timestep;
      PotentialEnergy pe(radius, data_recv, decomp_str);
      pe.aggregate_per_cell();
      auto const properties_map = pe.getAggregateMap();
      aggr_result_t aggregate{.agg = properties_map, .ts = timestep};
      world.send(0, TAG_RESULT, aggregate);
      world.recv(0, TAG_BREAK, stop);
      BOOST_LOG_TRIVIAL(info) << "[SLAVE: " << rank << "] got properties)";
    }
  }

  return EXIT_SUCCESS;
}
void submitPieces(const boost::mpi::communicator &world, Config &runningConf,
                  int t_len, std::deque<long> &ts,
                  decom_vec_storage_t decomp_str,
                  std::vector<aggr_result_t> &results,
                  const boost::filesystem::path &system_path,
                  boost::filesystem::path &cellstr_path, c_storage_index_t &db,
                  int &chunk_len, long &v_index,
                  std::vector<boost::mpi::request> &reqs_world, bool &stop,
                  int world_size, long &counter)
{
  while (!ts.empty()) {
    if (auto f_rank = test_any(reqs_world.begin(), reqs_world.end())) {
      auto dst_rank = f_rank.value().first.source();
      auto iterator = f_rank.value().second;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] send Job to Rank " << dst_rank << "\n";

      world.send(dst_rank, TAG_BREAK, stop);

      // Send the new job.
      long timestep{ts.front()};
      auto cols = db.select(
          columns(&ContactColumns::p1_id, &ContactColumns::p2_id,
                  &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
                  &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
                  &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
                  &ContactColumns::ct_force_z, &ContactColumns::cellstr,
                  &ContactColumns::ts),
          where(c(&ContactColumns::ts) == timestep));
      auto col_size = cols.size();
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending new job (" << timestep
                              << ") to SLAVE " << dst_rank;
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      counter++;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (static_cast<double>(counter) / t_len) * 100.0
          << "% done";
      world.send(dst_rank, TAG_SIZE, col_size);
      world.send(dst_rank, TAG_FILE, cols);
      ts.pop_front();
      *iterator = world.irecv(dst_rank, TAG_RESULT, results[v_index++]);
      if (v_index == chunk_len) {
        BOOST_LOG_TRIVIAL(info)
            << "[MASTER] Vector capacity exceeded. Write results.\n";
        wait_all(reqs_world.begin(), reqs_world.end());
        write_results(runningConf, decomp_str, results, system_path,
                      cellstr_path);
        v_index = 0;
        reqs_world.clear();
        results.clear();
        results.resize(chunk_len);
        for (dst_rank = 1; dst_rank < world_size; ++dst_rank) {
          world.send(dst_rank, TAG_BREAK, stop);
        }
        submitCompleteWorld(world, world_size, t_len, ts, db, reqs_world,
                            results, v_index, counter);
      }
    }
  }
}
void submitCompleteWorld(const boost::mpi::communicator &world, int world_size,
                         int t_len, std::deque<long> &ts, c_storage_index_t &db,
                         std::vector<boost::mpi::request> &reqs_world,
                         std::vector<aggr_result_t> &results, long &v_index,
                         long &counter)
{
  for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
    long timestep{ts.front()};
    auto cols = db.select(
        columns(&ContactColumns::p1_id, &ContactColumns::p2_id,
                &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
                &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
                &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
                &ContactColumns::ct_force_z, &ContactColumns::cellstr,
                &ContactColumns::ts),
        where(c(&ContactColumns::ts) == timestep));
    auto col_size = cols.size();
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job " << timestep
                            << " to SLAVE (first loop) " << dst_rank << "\n";
    BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
    counter++;
    BOOST_LOG_TRIVIAL(info)
    << "[MASTER] " << (static_cast<double>(counter) / t_len) * 100.0
    << "% done";
    world.send(dst_rank, TAG_SIZE, col_size);
    world.send(dst_rank, TAG_FILE, cols);
    // Post receive request for new jobs requests by slave [nonblocking]
    reqs_world.push_back(world.irecv(dst_rank, TAG_RESULT, results[v_index++]));
    ts.pop_front();
    if (ts.empty())
      break;
  }
}

void write_results(Config &runningConf, decom_vec_storage_t decomp_str,
                   std::vector<aggr_result_t> &results,
                   const boost::filesystem::path &system_path,
                   boost::filesystem::path &cellstr_path)
{
  BOOST_LOG_TRIVIAL(info) << "Writing system sum";
  std::ofstream f(system_path.string() + "/system.csv", std::ios::app);
  calc_system_sum_write_file(results, f, runningConf);
  BOOST_LOG_TRIVIAL(info) << "Finished writing system sum";

  BOOST_LOG_TRIVIAL(info) << "Writing individual particles";
  BOOST_LOG_TRIVIAL(info) << "Create cellstr map";
  BOOST_LOG_TRIVIAL(info) << "Intitalize cellfiles";
  BOOST_LOG_TRIVIAL(info) << "Intitalize cellfiles";
  write_cellstr_res(decomp_str, results, runningConf, cellstr_path);
  BOOST_LOG_TRIVIAL(info) << "Finished cellstr map";
}
void initialize_output_files(const Config &runningConf,
                             decom_vec_storage_t decomp_str,
                             const boost::filesystem::path &system_path,
                             const boost::filesystem::path &cellstr_path)
{
  std::ofstream f(system_path.string() + "/system.csv");
  f << "ts" << runningConf.sep << "penor" << runningConf.sep << "petan"
    << runningConf.sep << "ftan" << runningConf.sep << "fnor\n";
  f.close();
  BOOST_LOG_TRIVIAL(info) << "Intitalize cellfiles";
  for (auto &elem : decomp_str) {
    std::ofstream f{cellstr_path.string() + "/" + elem.cellstr + ".csv"};
    f << "ts" << runningConf.sep << "penor" << runningConf.sep << "petan"
      << runningConf.sep << "ftan" << runningConf.sep << "fnor\n";
  }
}
void write_cellstr_res(decom_vec_storage_t decomp_str,
                       std::vector<aggr_result_t> &results, Config &runningConf,
                       boost::filesystem::path &cellstr_path)
{
  for (auto &elem : decomp_str) {
    std::ofstream f(cellstr_path.string() + "/" + elem.cellstr + ".csv",
                    std::ios::app);
    for (auto &ts : results) {
      if (ts.ts == 0)
        break;
      f << ts.ts << runningConf.sep;
      f << ts.agg[elem.cellstr]["penor"] << runningConf.sep;
      f << ts.agg[elem.cellstr]["petan"] << runningConf.sep;
      f << ts.agg[elem.cellstr]["ftan"] << runningConf.sep;
      f << ts.agg[elem.cellstr]["fnor"] << "\n";
    }
  }
}
void calc_system_sum_write_file(std::vector<aggr_result_t> &results,
                                std::ofstream &f, Config &runningConf)
{
  for (auto &elem : results) {
    Accumulator<double> penor, petan, ftan, fnor;
    if (elem.ts == 0)
      break;
    for (auto &cell : elem.agg) {
      penor(cell.second["penor"]);
      petan(cell.second["petan"]);
      ftan(cell.second["ftan"]);
      fnor(cell.second["fnor"]);
    }
    f << elem.ts << runningConf.sep << penor.getAccVal() << runningConf.sep
      << petan.getAccVal() << runningConf.sep << ftan.getAccVal()
      << runningConf.sep << fnor.getAccVal() << "\n";
  }
}