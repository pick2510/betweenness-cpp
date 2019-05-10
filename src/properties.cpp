//
// Created by std on 10.05.19.
//

#include "properties.h"
#include "INIReader.h"
#include "data.h"
#include "grid.h"
#include "sqlite_orm.h"
#include "utils.h"
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
#include <map>
#include <serialize_tuple.h>

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
  std::vector<int> keys{};
  std::vector<double> vals{};
  std::vector<t_ts_pot_res> results{};
  gethostname(hostname, HOSTNAME_LEN);
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
    auto lookup_table = get_lookup_table(radius_map);

    keys = getKeys(radius_map);
    vals = getVals(radius_map);
  }

  broadcast(world, keys, MASTER);
  broadcast(world, vals, MASTER);
  broadcast(world, path, MASTER);

  if (rank == MASTER) {
    // sanity check of file list
    t_len = ts.size();
    BOOST_LOG_TRIVIAL(info) << "t_len = " << t_len;
    results.resize(t_len);
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
    std::vector<boost::mpi::request> reqs_world(world_size);

    for (int dst_rank = 1; dst_rank < world_size; ++dst_rank) {
      long timestep{ts.front()};
      auto db = indexContactStorage(path);
      db.sync_schema();
      auto cols = db.select(
          columns(&ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
                  &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
                  &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
                  &ContactColumns::ct_force_z, &ContactColumns::cellstr),
          where(c(&ContactColumns::ts) == timestep));
      auto col_size = cols.size();
      BOOST_LOG_TRIVIAL(info) << "[MASTER] Sending job " << timestep
                              << " to SLAVE (first loop) " << dst_rank << "\n";
      BOOST_LOG_TRIVIAL(info) << "[MASTER] v_index = " << v_index;
      BOOST_LOG_TRIVIAL(info)
          << "[MASTER] " << (v_index / t_len) * 100.0 << "% done";
      world.send(dst_rank, TAG_SIZE, col_size);
      world.send(dst_rank, TAG_FILE, cols);
      ts.pop_front();
      // Post receive request for new jobs requests by slave [nonblocking]
      reqs_world[dst_rank - 1] =
          world.irecv(dst_rank, TAG_RESULT, results[v_index++]);
    }
  }
  return EXIT_SUCCESS;
}