#include "graph_sqlite.h"

GraphSQLite::GraphSQLite(const std::map<std::string, int> &vertices_map,
                         const std::string &stor, long timestep)
    : v_map{vertices_map}, path{stor}, timestep{timestep}, graph{v_map.size()}
{
  BOOST_LOG_TRIVIAL(info) << timestep;
  BOOST_LOG_TRIVIAL(info) << "Initialized";
}

void GraphSQLite::calc()
{
  generate_graph();
  calculate_betweenness_centrality();
  calculate_accumulator();
}

void GraphSQLite::calculate_accumulator()
{

  for (auto &val : vals) {

    acc(val);
  }
  BOOST_LOG_TRIVIAL(info) << "Mean " << boost::accumulators::mean(acc);
}

void GraphSQLite::generate_graph()
{
  auto db = indexStorage(path);
  auto columns =
      db.get_all<ContactColumns>(where(c(&ContactColumns::ts) == timestep));
  for (const auto &elem : columns) {
    boost::add_edge(v_map.at(std::to_string(elem.p1_id)),
                    v_map.at(std::to_string(elem.p2_id)), graph);
  }
  BOOST_LOG_TRIVIAL(info) << "Graph finished!";
}
void GraphSQLite::calculate_betweenness_centrality()
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

std::map<int, double> GraphSQLite::get_centrality_map()
{
  return constructMap(keys, vals);
}

long GraphSQLite::get_timestep() { return timestep; }

double GraphSQLite::get_mean() { return boost::accumulators::mean(acc); }

Result GraphSQLite::get_result()
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
GraphSQLite::~GraphSQLite() {}