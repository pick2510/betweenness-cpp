#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/betweenness_centrality.hpp>

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Dump_Graph;
typedef boost::shared_array_property_map<double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type> Centrality_Map;


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
    std::map<std::string,int> &v_map;
    Dump_Graph graph;
    long timestep;
    void set_fpointer(int n);
  
  public:
    Graph(std::string &Path, std::map<std::string,int> &vertices_map);
    void generate_graph();
    void calculate_betweenness_centrality();
    ~Graph();
};

#endif