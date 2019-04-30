#include <regrid.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "INIReader.h"
#include "data.h"
#include "decomposition.h"
#include "dumpfile.h"
#include "natural_sort.hpp"
#include "sqlite_orm.h"
#include "utils.h"
#include <Eigen/Eigen>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>

using namespace sqlite_orm;

int main(int argc, char **argv)
{
  char hostname[HOSTNAME_LEN]{};
  std::string configPath{};
  INIReader reader;
  gethostname(hostname, HOSTNAME_LEN);

  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Gridding of DEM Output files";
  BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
  BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
  try {
    configPath = getConfigPath(argc, argv);
  }
  catch (std::invalid_argument e) {
    BOOST_LOG_TRIVIAL(error) << e.what();
    exit(EXIT_FAILURE);
  }
  try {
    reader = parseConfigFile(configPath);
  }
  catch (std::invalid_argument e) {
    BOOST_LOG_TRIVIAL(error) << e.what();
    exit(EXIT_FAILURE);
  }
  Config runningConf = getGridConfigObj(reader);
  LogConfig(runningConf);
  const Decomposition decomp(runningConf);
  const std::string path{runningConf.OutputPath + "/" +
                         runningConf.db_filename};
  std::map<int, double> radius_map;
  auto stor = inittsstorage(path);
  auto particles = indexStorage(path);
  particles.sync_schema();
  auto radius_storage = initRadstorage(path);
  for (auto &elem : radius_storage.iterate<radius>()) {
    radius_map[elem.id] = elem.rad;
  }
  auto ts_list = stor.get_all<ts_column>();
  for (auto const ts : ts_list) {
    auto selected_column = particles.get_all<ContactColumns>(
        where(c(&ContactColumns::ts) == ts.ts));
    particles.begin_transaction();
    for (auto &elem : selected_column) {
      auto coord = dumpfile::calc_contactpoint(elem, radius_map);
      auto cell = decomp.calc_cell_numeric(coord);
      elem.cell_x = cell.x;
      elem.cell_y = cell.y;
      elem.cell_z = cell.z;
      elem.cellstr = decomp.calc_cell(coord);
      particles.update(elem);
    }
    particles.commit();
    // particles.update(selected_column);
  }
  return EXIT_SUCCESS;
}
