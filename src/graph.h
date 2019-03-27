#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Dump_Graph; 
typedef boost::shared_array_property_map<double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type> Centrality_Map;

class Graph
{
    const int skip_lines = 9;
    const unsigned int particle_1 = 12;
    const unsigned int particle_2 = 13;

  private:
    std::ifstream file;
    std::vector<std::vector<std::string>> filecontent;
    Centrality_Map c_map;
    Dump_Graph graph;
    void advance_fpointer();
  
  public:
    Graph(std::string &Path);
    void generate_graph();
    void calculate_betweenness_centrality();
    ~Graph();
};

#endif