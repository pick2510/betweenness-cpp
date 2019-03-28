#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Dump_Graph;
typedef boost::vector_property_map<double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type> Centrality_Map;
typedef boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean>> Acc;

class Graph
{
  const unsigned int begin_line = 9;
  const unsigned int ts_line = 1;
  const unsigned int particle_1 = 12;
  const unsigned int particle_2 = 13;

private:
  std::ifstream file;
  std::vector<std::vector<std::string>> filecontent;
  Centrality_Map c_map;
  std::map<std::string, int> &v_map;
  std::map<int, double> betweenness_centrality;
  Dump_Graph graph;
  long timestep;
  Acc acc;
  void generate_graph();
  void calculate_betweenness_centrality();
  void calculate_accumulator();
  void set_fpointer(int n);

public:
  Graph(std::string &Path, std::map<std::string, int> &vertices_map);
  void calc();
  std::map<int, double> get_centrality_map();
  long get_timestep();
  double get_mean();
  ~Graph();
};

#endif