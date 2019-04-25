
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
#include "decomposition.h"
#include "dumpfile.h"
#include "grid.h"
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
  // MASTER CODE
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
  Decomposition decomp(runningConf);
  std::string chainpattern(runningConf.InputPath + "/postchain/*.chain");
  std::string xyzpattern(runningConf.InputPath + "/postxyz/*.tet");
  xyzpattern = trim(xyzpattern);
  chainpattern = trim(chainpattern);
  auto chain_file_list = glob_deq(chainpattern);
  auto radius_file_list = glob(xyzpattern);
  if (radius_file_list.size() < 1 || chain_file_list.size() < 1) {
    BOOST_LOG_TRIVIAL(error)
        << "Need at least one file in postchain and postxyz Dir.";
    exit(EXIT_FAILURE);
  }
  auto radius_file = read_radius_file(radius_file_list[0]);
  radius_file.pop_back();
  std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
  auto lookup_table = get_lookup_table(radius_file);
  SI::natural::sort(chain_file_list);
  auto chain_size = chain_file_list.size();

  using Storage = decltype(initStorage(""));
  Storage storage = initStorage(runningConf.OutputPath + "/DEM.db");
  storage.pragma.journal_mode(journal_mode::WAL);
  storage.pragma.synchronous(0);
  storage.sync_schema();
  std::atomic<long> index{0};
  double percent{0.0};
  std::vector<std::vector<ContactColumns>> chunk_res;
#if defined(_OPENMP)
  omp_lock_t mutex;
  omp_init_lock(&mutex);
#pragma omp parallel for ordered
#endif
  for (int i = 0; i < chain_size; i++) {
    dumpfile Dump(chain_file_list[i], lookup_table, decomp);
    Dump.parse_file();
#if defined(_OPENMP)
    omp_set_lock(&mutex);
#endif
    chunk_res.push_back(Dump.getData());
    index++;
    if (chunk_res.size() > 100) {
      BOOST_LOG_TRIVIAL(info) << "Starting Transaction";
      storage.transaction([&] {
        for (auto &column : chunk_res) {
          for (auto &res : column) {
            storage.insert(res);
          }
        }
        return true;
      });
      BOOST_LOG_TRIVIAL(info) << "Transaction Finished";
      chunk_res.clear();
    }
#if defined(_OPENMP)
    omp_unset_lock(&mutex);
#endif
    percent = (index / chain_size) * 100.0;
    BOOST_LOG_TRIVIAL(info)
        << percent << "% (" << index << " of " << chain_size << ") done";
  }
  storage.transaction([&] {
    for (auto &column : chunk_res) {
      for (auto &res : column) {
        storage.insert(res);
      }
    }
    return true;
  });
#if defined(_OPENMP)
  omp_destroy_lock(&mutex);
#endif
  BOOST_LOG_TRIVIAL(info) << "Finished insert, start with index";
  using Idx = decltype(indexStorage(""));
  Idx idx = indexStorage(runningConf.OutputPath + "/DEM.db");
  idx.pragma.journal_mode(journal_mode::WAL);
  idx.pragma.synchronous(0);
  idx.sync_schema();
  BOOST_LOG_TRIVIAL(info) << "Finished indexing";
  return EXIT_SUCCESS;
}

void LogConfig(Config &conf)
{
  BOOST_LOG_TRIVIAL(info) << "Using Input: " << conf.InputPath;
  BOOST_LOG_TRIVIAL(info) << "Using Output: " << conf.OutputPath;
  BOOST_LOG_TRIVIAL(info) << "x_cells: " << conf.x_cells;
  BOOST_LOG_TRIVIAL(info) << "y_cells: " << conf.y_cells;
  BOOST_LOG_TRIVIAL(info) << "z_cells: " << conf.z_cells;
  BOOST_LOG_TRIVIAL(info) << "Domain size: " << conf.domainsize_x << "x"
                          << conf.domainsize_y << "x" << conf.domainsize_z;
}