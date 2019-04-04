#include "graph.h"
#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/log/trivial.hpp>
#include <boost/accumulators/statistics.hpp>
#include "utils.h"
#include "data.h"

Graph::Graph(std::string &Path, const std::map<std::string, int> &vertices_map) :  v_map{vertices_map}, graph{vertices_map.size()}
{
    file.open(Path);
    set_fpointer(Graph::ts_line);
    file >> timestep;
    BOOST_LOG_TRIVIAL(info) << timestep;
    set_fpointer(Graph::begin_line);
    BOOST_LOG_TRIVIAL(info) << "Initialized";
}

void Graph::calc()
{
    generate_graph();
    calculate_betweenness_centrality();
    calculate_accumulator();
}

void Graph::calculate_accumulator()
{

    for (auto &kv : betweenness_centrality)
    {
        acc(kv.second);
    }
    BOOST_LOG_TRIVIAL(info) << "Mean " << boost::accumulators::mean(acc);
}

void Graph::set_fpointer(int n)
{
    goto_line(file, n);
}

void Graph::generate_graph()
{
    std::string line;
    while (std::getline(file, line))
    {
        std::vector<std::string> splitted_line;
        split_string(line, splitted_line);
        boost::add_edge(v_map.at(splitted_line[Graph::particle_1]), v_map.at(splitted_line[Graph::particle_2]), graph);
    }
    BOOST_LOG_TRIVIAL(info) << "Graph finished!";
}

void Graph::calculate_betweenness_centrality()
{
    c_map = Centrality_Map(boost::num_vertices(graph), boost::get(boost::vertex_index, graph));
    boost::brandes_betweenness_centrality(graph, c_map);
    BOOST_LOG_TRIVIAL(info) << "Calculated Betweenness Centrality";
    BGL_FORALL_VERTICES(vertex, graph, Dump_Graph)
    {
        auto val = c_map[vertex];
        betweenness_centrality[vertex] = val;
        v_betweeness.push_back(val);
    }
    std::sort(v_betweeness.begin(), v_betweeness.end());

}

std::map<int, double> Graph::get_centrality_map()
{
    return betweenness_centrality;
}

long Graph::get_timestep()
{
    return timestep;
}

double Graph::get_mean()
{
    return boost::accumulators::mean(acc);
}

Result Graph::get_result()
{
    Result res;
    res.b_centrality = betweenness_centrality;
    res.mean = accumulators::mean(acc);
    res.ts = timestep;
    res.q_090 = get_percentile_vector(v_betweeness, 0.90);
    res.q_099 = get_percentile_vector(v_betweeness, 0.99);
    res.kur = accumulators::kurtosis(acc);
    res.var = accumulators::variance(acc);
    res.skew = accumulators::skewness(acc);
    return res;
}

Graph::~Graph()
{
    file.close();
}
