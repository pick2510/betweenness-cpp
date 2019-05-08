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
  const std::string con_path{runningConf.OutputPath + "/" +
                             runningConf.contact_filename};
  const std::string part_path{runningConf.OutputPath + "/" +
                              runningConf.particle_filename};
  std::map<int, double> radius_map;
  auto ts_stor = initTSStorage(con_path);
  auto part_stor = ParticleIndexStorage(part_path);
  auto con_stor = indexContactStorage(con_path);
  auto radius_storage = initRadstorage(con_path);
  for (auto &elem : radius_storage.iterate<radius>()) {
    radius_map[elem.id] = elem.rad;
  }
  auto ts_list = ts_stor.get_all<ts_column>();
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp single
  {
#pragma omp task
    {
#endif
      for (const auto ts : ts_list) {
        BOOST_LOG_TRIVIAL(info) << "Contact ts: " << ts.ts;
        auto selected_column = con_stor.get_all<ContactColumns>(
            where(c(&ContactColumns::ts) == ts.ts));
        con_stor.begin_transaction();
        for (auto &elem : selected_column) {
          auto coord = dumpfile::calc_contactpoint(elem, radius_map);
          auto cell = decomp.calc_cell_numeric(coord);
          elem.cell_x = cell.x;
          elem.cell_y = cell.y;
          elem.cell_z = cell.z;
          elem.cellstr = decomp.calc_cell(coord);
          con_stor.update(elem);
        }
        con_stor.commit();
        // particles.update(selected_column);
      }
#if defined(_OPENMP)
    }
#pragma omp task
    {
#endif
      for (const auto ts : ts_list) {
        BOOST_LOG_TRIVIAL(info) << "Particle ts: " << ts.ts;
        auto selected_column = part_stor.get_all<ParticleColumns>(
            where(c(&ParticleColumns::ts) == ts.ts));
        part_stor.begin_transaction();
        for (auto &elem : selected_column) {
          auto coord = coordinate{.x = elem.p_x, .y = elem.p_y, .z = elem.p_z};
          auto cell = decomp.calc_cell_numeric(coord);
          elem.cell_x = cell.x;
          elem.cell_y = cell.y;
          elem.cell_z = cell.z;
          elem.cellstr = decomp.calc_cell(coord);
          part_stor.update(elem);
        }
        part_stor.commit();
        // particles.update(selected_column);
      }

#if defined(_OPENMP)
    }
  }
#endif
  return EXIT_SUCCESS;
}
