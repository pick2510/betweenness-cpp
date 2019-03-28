#include "graph.h"
#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include "utils.h"

Graph::Graph(std::string &Path, std::map<std::string,int> &vertices_map) : v_map(vertices_map), graph(vertices_map.size())
{
    file.open(Path);
    advance_fpointer();
    
}

void Graph::advance_fpointer(){
    std::string line;
    for (int i = 0; i < Graph::skip_lines; i++)
    {
        std::getline(file, line);
    }
}

void Graph::generate_graph(){
    std::string line;
    while (std::getline(file, line)){
        std::vector<std::string> splitted_line;
        split_string(line, splitted_line);
        boost::add_edge(v_map[splitted_line[Graph::particle_1]], v_map[splitted_line[Graph::particle_2]] , graph);

    }
    std::cout << "Graph finished!" << std::endl;

  
}

void Graph::calculate_betweenness_centrality(){
   //const double b_ij = 1 / ((v_map.size() - 1) * (v_map.size() - 2));
   auto c_map = Centrality_Map_Matrix(boost::num_vertices(graph), boost::get(boost::vertex_index, graph));
   boost::brandes_betweenness_centrality(graph, c_map);
   BGL_FORALL_VERTICES(vertex, graph, Dump_Graph){
       std::cout << vertex << ": " << c_map[vertex] << std::endl;
   }
}



Graph::~Graph()
{
    file.close();
}
