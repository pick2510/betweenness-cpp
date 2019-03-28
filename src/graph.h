#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <fstream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/betweenness_centrality.hpp>

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> Dump_Graph;
typedef boost::adjacency_matrix<boost::undirectedS> Dump_Graph_Matrix;
typedef boost::shared_array_property_map<double, boost::property_map<Dump_Graph, boost::vertex_index_t>::const_type> Centrality_Map;
typedef boost::shared_array_property_map<double, boost::property_map<Dump_Graph_Matrix, boost::vertex_index_t>::const_type> Centrality_Map_Matrix;


class Graph
{
    const int skip_lines = 9;
    const unsigned int particle_1 = 12;
    const unsigned int particle_2 = 13;

  private:
    std::ifstream file;
    std::vector<std::vector<std::string>> filecontent;
    Centrality_Map c_map;
    std::map<std::string,int> &v_map;
    Dump_Graph_Matrix graph;
    void advance_fpointer();
  
  public:
    Graph(std::string &Path, std::map<std::string,int> &vertices_map);
    void generate_graph();
    void calculate_betweenness_centrality();
    ~Graph();
};

#endif