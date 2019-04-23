
#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

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

  using Storage = decltype(initStorage(""));
  Storage storage = initStorage("");
  storage.sync_schema();
  return EXIT_SUCCESS;
}
