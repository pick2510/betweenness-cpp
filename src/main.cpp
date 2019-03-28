#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "main.h"
#include "utils.h"
#include "natural_sort.hpp"
#include "graph.h"

int main(int argc, char **argv)
{
  std::cout << "****************************************\n";
  std::cout << "Node betweenness centrality\n";
  std::cout << "Kintali (2008) arXiv:0809.1906v2 [cs.DS] \n";
  std::cout << "(generalization to directed and disconnected graphs)\n";
  std::cout << "C++ Version\n";
  std::cout << "by Dominik Strebel (2019)\n";
  std::cout << "****************************************\n";
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
  std::string chainpattern(runningConf.InputPath + "/*.chain");
  std::string xyzpattern(runningConf.InputPath + "/*.tet");
  xyzpattern = trim(xyzpattern);
  chainpattern = trim(chainpattern);
  std::vector<std::string> chain_file_list = glob(chainpattern);
  std::vector<std::string> radius_file_list = glob(xyzpattern);
  if (radius_file_list.empty())
  {
    std::cout << "Need at least one file with radii *.tet\n";
    exit(EXIT_FAILURE);
  }
  SI::natural::sort(chain_file_list);
  auto radius_file = read_file(radius_file_list[0]);
  radius_file.pop_back();
  std::sort(radius_file.begin(), radius_file.end(), cmp_radii);
  auto lookup_table = get_lookup_table(radius_file);
  auto vertice_map = get_vertice_map(radius_file);
  Graph mygraph(chain_file_list[0], vertice_map);
  mygraph.calc();
  auto mean = mygraph.get_mean();
}
