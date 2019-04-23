
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
#include <unistd.h>

#include "sqlite_orm.h"
#include  "grid.h"
#include "natural_sort.hpp"
#include "utils.h"
#include <Eigen/Eigen>
#include <boost/log/trivial.hpp>



int main(int argc, char **argv){
    char hostname[HOSTNAME_LEN]{};
    Config runningConf{};
    gethostname(hostname, HOSTNAME_LEN);
     // MASTER CODE
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Node betweenness centrality";
    BOOST_LOG_TRIVIAL(info) << "using BOOST Graph Library";
    BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
    BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
    BOOST_LOG_TRIVIAL(info) << "****************************************";
    BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
    try {
      runningConf = getCL(argc, argv);
    }
    catch (std::invalid_argument e) {
      BOOST_LOG_TRIVIAL(info) << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }
}