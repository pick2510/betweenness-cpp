
#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

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
  Config runningConf{};
  gethostname(hostname, HOSTNAME_LEN);
  // MASTER CODE
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Gridding of DEM Output files";
  BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
  BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
  try {
    runningConf = getCL(argc, argv);
  }
  catch (std::invalid_argument e) {
    BOOST_LOG_TRIVIAL(error) << e.what();
    exit(EXIT_FAILURE);
  }
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
  SI::natural::sort(chain_file_list);
  auto chain_size = chain_file_list.size();
  using Storage = decltype(initStorage(""));
  Storage storage = initStorage(runningConf.OutputPath + "/" + "DEM.db");
  storage.pragma.journal_mode(journal_mode::WAL);
  storage.pragma.synchronous(0);
  storage.sync_schema();
  std::atomic<long> index{0};
  double percent{0.0};
  std::vector<std::vector<ContactColumns>> chunk_res;
  omp_lock_t mutex;
  omp_init_lock(&mutex);
#pragma omp parallel for ordered
  for (int i = 0; i < chain_size; i++) {
    dumpfile Dump(chain_file_list[i]);
    Dump.parse_file();
    omp_set_lock(&mutex);
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
    omp_unset_lock(&mutex);
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
  omp_destroy_lock(&mutex);
  BOOST_LOG_TRIVIAL(info) << "Finished insert, start with index";
  auto idx = make_storage(runningConf.OutputPath + "/" + "DEM.db", make_index("idx_ts_cellstr", &ContactColumns::ts, &ContactColumns::cellstr));
  idx.sync_schema();
  BOOST_LOG_TRIVIAL(info) << "Finished indexing";
  return EXIT_SUCCESS;
}
