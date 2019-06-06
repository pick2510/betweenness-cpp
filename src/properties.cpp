//
// Created by std on 10.05.19.
//
#include "properties.h"
#include "INIReader.h"
#include "PropertyCalculator.h"
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
  std::string contact_path;
  std::string particle_path;
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
    BOOST_LOG_TRIVIAL(info) << "Property calculator";
    BOOST_LOG_TRIVIAL(info) << "for DEM Output, SQLite Version";
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

    contact_path =
        std::string{runningConf.InputPath + "/" + runningConf.contact_filename};

    particle_path = std::string{runningConf.InputPath + "/" +
                                runningConf.particle_filename};
    cellstr_path = boost::filesystem::path(runningConf.OutputPath + "/cells");
    system_path = boost::filesystem::path(runningConf.OutputPath + "/system");
    check_path(cellstr_path);
    check_path(system_path);
    auto stor = initTSStorage(contact_path);
    auto particles = indexContactStorage(contact_path);
    auto radius_storage = initRadstorage(contact_path);
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
    auto cell_db = initDecompstorage(contact_path);
    decomp_str = cell_db.get_all<decomp_table>();
    keys = getKeys(radius_map);
    vals = getVals(radius_map);
  }

  broadcast(world, keys, MASTER);
  broadcast(world, decomp_str, MASTER);
  broadcast(world, vals, MASTER);
  auto c_db = indexContactStorage(contact_path);
  auto p_db = ParticleIndexStorage(particle_path);
  if (rank == MASTER) {
    auto c_db = indexContactStorage(contact_path);
    auto p_db = ParticleIndexStorage(particle_path);
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

    submitCompleteWorld(world, world_size, t_len, ts, c_db, reqs_world, results,
                        v_index, counter, p_db);
    bool stop = false;

    if (!ts.empty()) {
      submitPieces(world, runningConf, t_len, ts, decomp_str, results,
                   system_path, cellstr_path, c_db, chunk_len, v_index,
                   reqs_world, stop, world_size, counter, p_db);
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sent all jobs.\n";
      wait_all(reqs_world.begin(), reqs_world.end());
    }
    else {
      wait_all(reqs_world.begin(), reqs_world.end());
    }
    initialize_output_files(runningConf, decomp_str, system_path, cellstr_path,
                            results[0]);
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
      tuple_contact_storage_t data_contact_recv{};
      tuple_particle_storage_t data_particle_recv{};
      size_t contact_col_size, particle_col_size;
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] (" << hostname << ") Intialized JOB";
      world.recv(0, TAG_SIZE, contact_col_size);
      data_contact_recv.resize(contact_col_size);
      world.recv(0, TAG_FILE, data_contact_recv);
      world.recv(0, TAG_SIZE, particle_col_size);
      data_particle_recv.resize(particle_col_size);
      world.recv(0, TAG_P_COL, data_particle_recv);
      auto timestep = std::get<14>(*data_contact_recv.begin());
      BOOST_LOG_TRIVIAL(info)
          << "[SLAVE: " << rank << "] recevied ts: " << timestep;
      PropertyCalculator pe(radius, data_contact_recv, data_particle_recv,
                            decomp_str);
      pe.aggregate_per_cell();
      auto const properties_map = pe.getAggregateMap();
      auto const global_map = pe.getGlobalMap();
      aggr_result_t aggregate{
          .agg = properties_map, .global = global_map, .ts = timestep};
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
                  boost::filesystem::path &cellstr_path,
                  c_storage_index_t &c_db, int &chunk_len, long &v_index,
                  std::vector<boost::mpi::request> &reqs_world, bool &stop,
                  int world_size, long &counter, p_storage_index_t &p_db)
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
      auto contact_cols = c_db.select(
          columns(&ContactColumns::p1_id, &ContactColumns::p2_id,
                  &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
                  &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
                  &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
                  &ContactColumns::ct_force_z, &ContactColumns::c_force_x,
                  &ContactColumns::c_force_y, &ContactColumns::c_force_z,
                  &ContactColumns::sliding_contact, &ContactColumns::cellstr,
                  &ContactColumns::ts),
          where(c(&ContactColumns::ts) == timestep));
      auto particle_cols = p_db.select(
          columns(&ParticleColumns::p_id, &ParticleColumns::p_vx,
                  &ParticleColumns::p_vy, &ParticleColumns::p_vz,
                  &ParticleColumns::p_coord, &ParticleColumns::p_disp_x,
                  &ParticleColumns::p_disp_y, &ParticleColumns::p_disp_z,
                  &ParticleColumns::p_ke_rot, &ParticleColumns::p_ke_tra,
                  &ParticleColumns::p_omegax, &ParticleColumns::p_omegay,
                  &ParticleColumns::p_omegaz, &ParticleColumns::cellstr,
                  &ParticleColumns::ts),
          where(c(&ParticleColumns::ts) == timestep));
      auto particle_col_size = particle_cols.size();
      auto contact_col_size = contact_cols.size();
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending new job (" << timestep
                              << ") to SLAVE " << dst_rank;
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      counter++;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (static_cast<double>(counter) / t_len) * 100.0
          << "% done";
      world.send(dst_rank, TAG_SIZE, contact_col_size);
      world.send(dst_rank, TAG_FILE, contact_cols);
      world.send(dst_rank, TAG_SIZE, particle_col_size);
      world.send(dst_rank, TAG_P_COL, particle_cols);
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
        submitCompleteWorld(world, world_size, t_len, ts, c_db, reqs_world,
                            results, v_index, counter, p_db);
      }
    }
  }
}
void submitCompleteWorld(const boost::mpi::communicator &world, int world_size,
                         int t_len, std::deque<long> &ts,
                         c_storage_index_t &c_db,
                         std::vector<boost::mpi::request> &reqs_world,
                         std::vector<aggr_result_t> &results, long &v_index,
                         long &counter, p_storage_index_t &p_db)
{
  for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
    long timestep{ts.front()};
    auto contact_cols = c_db.select(
        columns(&ContactColumns::p1_id, &ContactColumns::p2_id,
                &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
                &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
                &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
                &ContactColumns::ct_force_z, &ContactColumns::c_force_x,
                &ContactColumns::c_force_y, &ContactColumns::c_force_z,
                &ContactColumns::sliding_contact, &ContactColumns::cellstr,
                &ContactColumns::ts),
        where(c(&ContactColumns::ts) == timestep));
    auto particle_cols = p_db.select(
        columns(&ParticleColumns::p_id, &ParticleColumns::p_vx,
                &ParticleColumns::p_vy, &ParticleColumns::p_vz,
                &ParticleColumns::p_coord, &ParticleColumns::p_disp_x,
                &ParticleColumns::p_disp_y, &ParticleColumns::p_disp_z,
                &ParticleColumns::p_ke_rot, &ParticleColumns::p_ke_tra,
                &ParticleColumns::p_omegax, &ParticleColumns::p_omegay,
                &ParticleColumns::p_omegaz, &ParticleColumns::cellstr,
                &ParticleColumns::ts),
        where(c(&ParticleColumns::ts) == timestep));
    auto col_size = contact_cols.size();
    auto particle_col_size = particle_cols.size();
    BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job " << timestep
                            << " to SLAVE (first loop) " << dst_rank << "\n";
    BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
    counter++;
    BOOST_LOG_TRIVIAL(info)
        << "[MASTER] " << (static_cast<double>(counter) / t_len) * 100.0
        << "% done";
    world.send(dst_rank, TAG_SIZE, col_size);
    world.send(dst_rank, TAG_FILE, contact_cols);
    world.send(dst_rank, TAG_SIZE, particle_col_size);
    world.send(dst_rank, TAG_P_COL, particle_cols);
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
  std::ofstream f(system_path.string() + "/global.csv", std::ios::app);
  write_global_state(results, f, runningConf);
  BOOST_LOG_TRIVIAL(info) << "Finished writing system sum";

  BOOST_LOG_TRIVIAL(info) << "Writing individual particles";
  BOOST_LOG_TRIVIAL(info) << "Create cellstr map";
  BOOST_LOG_TRIVIAL(info) << "Intitalize cellfiles";
  write_cellstr_res(decomp_str, results, runningConf, cellstr_path);
  BOOST_LOG_TRIVIAL(info) << "Finished cellstr map";
}
void initialize_output_files(const Config &runningConf,
                             decom_vec_storage_t decomp_str,
                             const boost::filesystem::path &system_path,
                             const boost::filesystem::path &cellstr_path,
                             aggr_result_t &res)
{
  std::ofstream f(system_path.string() + "/global.csv");
  f << "ts" << runningConf.sep;
  auto cols_size = sizeof(f_header) / sizeof(f_header[0]);
  for (int i = 0; i < cols_size; i++) {
    f << f_header[i] << runningConf.sep;
  }
  f << "\n";
  BOOST_LOG_TRIVIAL(info) << "Intitalize cellfiles";
  for (auto &elem : decomp_str) {
    std::ofstream f{cellstr_path.string() + "/" + elem.cellstr + ".csv"};
    f << "ts" << runningConf.sep;
    for (int i = 0; i < cols_size; i++) {
      f << f_header[i] << runningConf.sep;
    }
    f << "\n";
  }
}
void write_cellstr_res(decom_vec_storage_t decomp_str,
                       std::vector<aggr_result_t> &results, Config &runningConf,
                       boost::filesystem::path &cellstr_path)
{
  for (auto &elem : decomp_str) {
    std::ofstream f(cellstr_path.string() + "/" + elem.cellstr + ".csv",
                    std::ios::app);
    auto cols_size = sizeof(f_header) / sizeof(f_header[0]);
    for (auto &ts : results) {
      if (ts.ts == 0)
        break;
      f << ts.ts << runningConf.sep;
      for (int i = 0; i < cols_size; i++) {
        // for (auto &cell : ts.agg[elem.cellstr]) {
        f << ts.agg[elem.cellstr][f_header[i]] << runningConf.sep;
      }
      f << "\n";
    }
  }
}
void write_global_state(std::vector<aggr_result_t> &results, std::ofstream &f,
                        Config &runningConf)
{
  auto cols_size = sizeof(f_header) / sizeof(f_header[0]);
  for (auto &ts : results) {
    if (ts.ts == 0)
      break;
    f << ts.ts << runningConf.sep;
    for (int i = 0; i < cols_size; i++) {
      f << ts.global[f_header[i]] << runningConf.sep;
    }
    f << "\n";
  }
}