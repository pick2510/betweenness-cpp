#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <fstream>
#include <vector>
#include "data.h"
#include <boost/array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace tags = boost::accumulators::tag;
namespace accumulators = boost::accumulators;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Dump_Graph;
typedef boost::vector_property_map<double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type> Centrality_Map;
typedef accumulators::accumulator_set<double, accumulators::features<tags::mean,tags::variance, tags::skewness, tags::kurtosis>> Acc;

class Graph
{
  static constexpr int begin_line = 9;
  static constexpr int ts_line = 1;
  static constexpr int particle_1 = 12;
  static constexpr int particle_2 = 13;

private:
  std::ifstream file;
  std::vector<std::vector<std::string>> filecontent;
  Centrality_Map c_map;
  const std::map<std::string, int> &v_map;
  //std::map<int, double> betweenness_centrality;
  std::vector<int> keys;
  std::vector<double> vals;
  std::vector<double> v_betweeness;
  Dump_Graph graph;
  long timestep;
  Acc acc;
  void generate_graph();
  void calculate_betweenness_centrality();
  void calculate_accumulator();
  void set_fpointer(int n);

public:
  Graph(std::string &Path, const std::map<std::string, int> &vertices_map);
  Result get_result();
  void calc();
  std::map<int, double> get_centrality_map();
  long get_timestep();
  double get_mean();
  ~Graph();
};


template <typename K,typename V>
V get_percentile_map(std::map<K,V> map, double percentile){
    std::vector<V> valVec;
    for (auto &kv : map){
        valVec.push_back(kv.second);
    }
    std::sort(valVec.begin(), valVec.end());
    int n = std::ceil(percentile * valVec.size()) - 1;
    return valVec[n];
}


template <typename T>
T get_percentile_vector(std::vector<T> &valVec, double percentile){
  int n = std::round(percentile * valVec.size()) - 1;
  return valVec[n];
}

#endif