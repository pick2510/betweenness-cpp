#include "graph.h"
#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/log/trivial.hpp>
#include "utils.h"
#include "data.h"


Graph::Graph(std::string &Path, std::map<std::string,int> &vertices_map) : v_map(vertices_map), graph(vertices_map.size())
{
    file.open(Path);
    set_fpointer(Graph::ts_line);
    file >> timestep;
    BOOST_LOG_TRIVIAL(info) << timestep;
    set_fpointer(Graph::begin_line);
    
}


void Graph::set_fpointer(int n){
    goto_line(file,n);
}


void Graph::generate_graph(){
    std::string line;
    while (std::getline(file, line)){
        std::vector<std::string> splitted_line;
        split_string(line, splitted_line);
        boost::add_edge(v_map[splitted_line[Graph::particle_1]], v_map[splitted_line[Graph::particle_2]] , graph);

    }
    BOOST_LOG_TRIVIAL(info) << "Graph finished!";

  
}

void Graph::calculate_betweenness_centrality(){
   //const double b_ij = 1 / ((v_map.size() - 1) * (v_map.size() - 2));
   auto c_map = Centrality_Map(boost::num_vertices(graph), boost::get(boost::vertex_index, graph));
   boost::brandes_betweenness_centrality(graph, c_map);
   BGL_FORALL_VERTICES(vertex, graph, Dump_Graph){
       std::cout << vertex << ": " << c_map[vertex] << std::endl;
   }
}



Graph::~Graph()
{
    file.close();
}
