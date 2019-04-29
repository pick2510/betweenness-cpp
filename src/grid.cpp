
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
  std::vector<double> lookup_table;
  std::map<int, double> radius_map;
  get_lookup_table(radius_file, lookup_table, radius_map);
  SI::natural::sort(chain_file_list);
  auto chain_size = chain_file_list.size();
  auto storage =
      initStorage(runningConf.OutputPath + "/" + runningConf.db_filename);
  storage.pragma.journal_mode(journal_mode::WAL);
  storage.pragma.synchronous(0);
  storage.sync_schema();
  std::atomic<long> index{0};
  double percent{0.0};
  std::vector<std::vector<ContactColumns>> chunk_res;
  std::vector<long> ts_vec;
#if defined(_OPENMP)
  omp_lock_t mutex;
  omp_init_lock(&mutex);
#pragma omp parallel for ordered
#endif
  for (int i = 0; i < chain_size; i++) {
    dumpfile Dump(chain_file_list[i], radius_map, decomp);
    Dump.parse_file();
#if defined(_OPENMP)
    omp_set_lock(&mutex);
#endif
    ts_vec.push_back(Dump.gettimestep());
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
    for (const auto &column : chunk_res) {
      for (const auto &res : column) {
        storage.insert(res);
      }
    }
    return true;
  });
#if defined(_OPENMP)
  omp_destroy_lock(&mutex);
#endif
  BOOST_LOG_TRIVIAL(info) << "Finished insert, start with index";
  auto idx =
      indexStorage(runningConf.OutputPath + "/" + runningConf.db_filename);
  idx.pragma.journal_mode(journal_mode::WAL);
  idx.pragma.synchronous(0);
  idx.sync_schema();
  BOOST_LOG_TRIVIAL(info) << "Finished indexing";
  std::sort(ts_vec.begin(), ts_vec.end());
  BOOST_LOG_TRIVIAL(info) << "Starting insert TS Table";
  auto ts_table =
      inittsstorage(runningConf.OutputPath + "/" + runningConf.db_filename);
  ts_table.sync_schema();
  ts_table.transaction([&] {
    for (const auto &ts : ts_vec) {
      ts_table.insert(ts_column{.ts = ts});
    }
    return true;
  });
  BOOST_LOG_TRIVIAL(info) << "Finished TS Table";
  BOOST_LOG_TRIVIAL(info) << "Starting insert Radius Table";
  auto rad_table =
      initRadstorage(runningConf.OutputPath + "/" + runningConf.db_filename);
  rad_table.sync_schema();
  rad_table.transaction([&] {
    for (const auto &r : radius_map) {
      rad_table.insert(radius{.id = r.first, .rad = r.second});
    }
    return true;
  });
  BOOST_LOG_TRIVIAL(info) << "Finished insert TS Table";
  return EXIT_SUCCESS;
}
