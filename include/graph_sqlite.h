#ifndef GRAPH_SQLITE_H
#define GRAPH_SQLITE_H

#include "data.h"
#include "grid.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/log/trivial.hpp>
#include <fstream>
#include <string>
#include <vector>

using Storage = decltype(indexContactStorage(""));
using namespace sqlite_orm;
namespace tags = boost::accumulators::tag;
namespace accumulators = boost::accumulators;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
    Dump_Graph;
typedef boost::vector_property_map<
    double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type>
    Centrality_Map;
typedef accumulators::accumulator_set<
    double, accumulators::features<tags::mean, tags::variance, tags::skewness,
                                   tags::kurtosis>>
    Acc;

class GraphSQLite {
  static constexpr int begin_line = 9;
  static constexpr int ts_line = 1;

private:
  Centrality_Map c_map;
  std::map<std::string, int> &v_map;
  long timestep;
  std::vector<ContactColumns> contact_cols;
  // Storage db;
  // std::map<int, double> betweenness_centrality;
  std::vector<int> keys;
  std::vector<double> vals;
  std::vector<double> v_betweeness;
  std::vector<ParticleColumns> &part_cols;
  std::multimap<std::string, double> cell_agg;
  Dump_Graph graph;
  const std::vector<decomp_table> &decomp;
  std::map<int, std::string> inv_map, part_cell_map;
  Acc acc_total_domain;
  void generate_graph();
  void calculate_betweenness_centrality();
  void calculate_accumulators();

public:
  GraphSQLite(std::map<std::string, int> &vertices_map, long ts,
              std::vector<ContactColumns> &contact_cols,
              const std::vector<decomp_table> &decomp,
              std::vector<ParticleColumns> &part_cols);
  Betweenness_Result get_result();
  void calc();
  std::map<int, double> get_centrality_map();
  long get_timestep();
  double get_mean();
  ~GraphSQLite();
};

template <typename K, typename V>
V get_percentile_map(std::map<K, V> map, double percentile)
{
  std::vector<V> valVec;
  for (auto &kv : map) {
    valVec.push_back(kv.second);
  }
  std::sort(valVec.begin(), valVec.end());
  int n = std::ceil(percentile * valVec.size()) - 1;
  return valVec[n];
}

template <typename T>
T get_percentile_vector(std::vector<T> &valVec, double percentile)
{
  int n = std::round(percentile * valVec.size()) - 1;
  return valVec[n];
}

#endif