#include "graph_sqlite.h"
#include "data.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <string>

using namespace sqlite_orm;

template <class T>
GraphSQLite<T>::GraphSQLite(const std::map<std::string, int> &vertices_map,
                            std::string &storage, long timestep)
    : v_map{vertices_map}, graph{vertices_map.size()}, timestep{timestep}
{
  std::string path{storage};
  T db = initStorage(path);
  BOOST_LOG_TRIVIAL(info) << timestep;
  BOOST_LOG_TRIVIAL(info) << "Initialized";
}

template <class T> void GraphSQLite<T>::calc()
{
  generate_graph();
  calculate_betweenness_centrality();
  calculate_accumulator();
}

template <class T> void GraphSQLite<T>::calculate_accumulator()
{

  for (auto &val : vals) {
    acc(val);
  }
  BOOST_LOG_TRIVIAL(info) << "Mean " << boost::accumulators::mean(acc);
}

template <class T> void GraphSQLite<T>::generate_graph()
{
  auto columns = db.template get_all<ContactColumns>(
      where(c(&ContactColumns::ts) == timestep));
  for (const auto &elem : columns) {
    boost::add_edge(v_map.at(elem.p1_id), v_map.at(elem.p2_id), graph);
  }
  BOOST_LOG_TRIVIAL(info) << "Graph finished!";
}

template <class T> void GraphSQLite<T>::calculate_betweenness_centrality()
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

template <class T> std::map<int, double> GraphSQLite<T>::get_centrality_map()
{
  return constructMap(keys, vals);
}

template <class T> long GraphSQLite<T>::get_timestep() { return timestep; }

template <class T> double GraphSQLite<T>::get_mean()
{
  return boost::accumulators::mean(acc);
}

template <class T> Result GraphSQLite<T>::get_result()
{
  Result res;
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

template <class T> GraphSQLite<T>::~GraphSQLite() {}
