#include "graph.h"
#include "data.h"
#include "utils.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <string>

Graph::Graph(std::string &Path, const std::map<std::string, int> &vertices_map)
    : v_map{vertices_map}, graph{vertices_map.size()}
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

  for (auto &val : vals) {
    acc(val);
  }
  BOOST_LOG_TRIVIAL(info) << "Mean " << boost::accumulators::mean(acc);
}

void Graph::set_fpointer(int n) { goto_line(file, n); }

void Graph::generate_graph()
{
  int error_count = 0;
  std::string line;
  while (std::getline(file, line)) {
    std::vector<std::string> splitted_line;
    split_string(line, splitted_line);
    try {
      boost::add_edge(v_map.at(splitted_line[Graph::particle_1]),
                      v_map.at(splitted_line[Graph::particle_2]), graph);
    }
    catch (const std::exception &e) {
      BOOST_LOG_TRIVIAL(error)
          << "GRAPH ERROR: Edge between: " << splitted_line[Graph::particle_1]
          << " and " << splitted_line[Graph::particle_2]
          << ". Didn't find in radius file";
      error_count++;
    }
  }
  BOOST_LOG_TRIVIAL(info) << "Graph finished!";
  if (error_count != 0) {
    BOOST_LOG_TRIVIAL(error)
        << "ERROR: Timestep: " << timestep << " Count: " << error_count << "\n";
  }
}

void Graph::calculate_betweenness_centrality()
{
  c_map = Centrality_Map(boost::num_vertices(graph),
                         boost::get(boost::vertex_index, graph));
  boost::brandes_betweenness_centrality(graph, c_map);
  BOOST_LOG_TRIVIAL(info) << "Calculated Betweenness Centrality";
  BGL_FORALL_VERTICES(vertex, graph, Dump_Graph)
  {
    auto val = c_map[vertex];
    keys.push_back(vertex);
    vals.push_back(val);
  }
  v_betweeness.resize(vals.size());
  std::copy(vals.begin(), vals.end(), v_betweeness.begin());
  std::sort(v_betweeness.begin(), v_betweeness.end());
}

std::map<int, double> Graph::get_centrality_map()
{
  return constructMap(keys, vals);
}

long Graph::get_timestep() { return timestep; }

double Graph::get_mean() { return boost::accumulators::mean(acc); }

Betweenness_Result Graph::get_result()
{
  Betweenness_Result res;
  res.keys = keys;
  res.vals = vals;
  res.mean = accumulators::mean(acc);
  res.ts = timestep;
  res.q_090 = get_percentile_vector(v_betweeness, 0.90);
  res.q_099 = get_percentile_vector(v_betweeness, 0.99);
  res.kur = accumulators::kurtosis(acc);
  res.var = accumulators::variance(acc);
  res.skew = accumulators::skewness(acc);
  return res;
}

Graph::~Graph() { file.close(); }
