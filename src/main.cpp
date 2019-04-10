#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <atomic>

#include "main.h"
#include "utils.h"
#include "natural_sort.hpp"
#include "graph.h"
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include "omp.h"

int main(int argc, char **argv)
{
  //boost::log::add_console_log(std::cout, boost::log::keywords::format = "[%TimeStamp%] [%Severity%] %Message%");
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Node betweenness centrality";
  BOOST_LOG_TRIVIAL(info) << "using BOOST Graph Library";
  BOOST_LOG_TRIVIAL(info) << "Calculates different properties";
  BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Using [" << omp_get_max_threads() << "] threads";
  Config runningConf;
  try
  {
    runningConf = getCL(argc, argv);
  }
  catch (std::invalid_argument e)
  {
    std::cout << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

// define pattern fro globbing chain files, trim strings, and glob

  std::string chainpattern(runningConf.InputPath + "/*.chain");
  std::string xyzpattern(runningConf.InputPath + "/*.tet");
  xyzpattern = trim(xyzpattern);
  chainpattern = trim(chainpattern);
  std::vector<std::string> chain_file_list = glob(chainpattern);
  std::vector<std::string> radius_file_list = glob(xyzpattern);

// sanity check of file list 

  if (radius_file_list.empty())
  {
    std::cout << "Need at least one file with radii *.tet\n";
    exit(EXIT_FAILURE);
  }
  SI::natural::sort(chain_file_list);

// Read one radius file, create lookup table for continous particle ids.
// Create 
  auto radius_file = read_file(radius_file_list[0]);
  radius_file.pop_back();
  std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
  auto lookup_table = get_lookup_table(radius_file);
  auto vertice_map = get_vertice_map(radius_file);
  auto inv_vertice_map = inverse_map(vertice_map);
  auto t_len = chain_file_list.size();
  
  
  std::vector<Result> results;
  omp_lock_t mutex;
  omp_init_lock(&mutex);
  std::atomic<long> index {0};
  double percent {0.0};
#pragma omp parallel for
  for (std::size_t i = 0; i < chain_file_list.size(); ++i)
  {
    Graph mygraph(chain_file_list[i], vertice_map);
    mygraph.calc();
    omp_set_lock(&mutex);
    results.push_back(mygraph.get_result());
    index++;
    omp_unset_lock(&mutex);
    percent = (index / t_len) * 100;
    BOOST_LOG_TRIVIAL(info) << percent << "% (" << index << " of " << t_len << ") done";
  }
  omp_destroy_lock(&mutex);
  std::sort(results.begin(), results.end(), cmp_ts);
  std::ofstream ts_mean_file(runningConf.OutputPath + "/properties.csv");
  write_ts_header(ts_mean_file, runningConf);
  ts_mean_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  for (auto &v : results)
  {
    ts_mean_file << std::to_string(v.ts) << runningConf.sep << std::to_string(v.mean)
                 << runningConf.sep << std::to_string(v.var) << runningConf.sep << std::to_string(std::sqrt(v.var))
                 << runningConf.sep << std::to_string(v.skew) << runningConf.sep << std::to_string(v.kur)
                 << runningConf.sep << std::to_string(v.q_090) << runningConf.sep << std::to_string(v.q_099)
                 << "\n";
    std::ofstream ts_file(runningConf.OutputPath + "/centrality_" + std::to_string(v.ts) + ".csv");
    write_cent_header(ts_file, runningConf);
    ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (auto &kv : v.b_centrality)
    {
      ts_file << inv_vertice_map[kv.first] << runningConf.sep << kv.second << "\n";
    }
    ts_file.close();
  }
  ts_mean_file.flush();
  ts_mean_file.close();
}
