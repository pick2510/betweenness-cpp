#include "graph_sqlite.h"

GraphSQLite::GraphSQLite(std::map<std::string, int> &vertices_map,
                         long timestep,
                         std::vector<ContactColumns> &contact_cols,
                         const std::vector<decomp_table> &decomp,
                         std::vector<ParticleColumns> &part_cols)
    : v_map{vertices_map}, timestep{timestep}, graph{v_map.size()},
      contact_cols{contact_cols}, decomp{decomp}, part_cols{part_cols}
{

  inv_map = inverse_map(v_map);
  for (auto &elem : part_cols) {
    part_cell_map[elem.p_id] = elem.cellstr;
  }
  BOOST_LOG_TRIVIAL(info) << timestep;
  BOOST_LOG_TRIVIAL(info) << "Initialized";
}

void GraphSQLite::calc()
{
  generate_graph();
  calculate_betweenness_centrality();
  calculate_accumulators();
}

void GraphSQLite::calculate_accumulators()
{

  for (int i = 0; i < vals.size(); i++) {

    acc_total_domain(vals[i]);
    cell_agg.insert(std::make_pair(part_cell_map.at(keys[i]), vals[i]));
  }
  BOOST_LOG_TRIVIAL(info) << "Mean "
                          << boost::accumulators::mean(acc_total_domain);
}

void GraphSQLite::generate_graph()
{
  // auto db = indexContactStorage(path);
  // auto columns =
  //   db.get_all<ContactColumns>(where(c(&ContactColumns::ts) == timestep));
  for (const auto &elem : contact_cols) {
    try {
      boost::add_edge(v_map.at(std::to_string(elem.p1_id)),
                      v_map.at(std::to_string(elem.p2_id)), graph);
    }
    catch (std::out_of_range &exc) {
      continue;
    }
    /*BOOST_LOG_TRIVIAL(info)
        << "GRAPH: p1_id-> " << elem.p1_id << ": "
        << v_map.at(std::to_string(elem.p1_id)) << " p2_id: " << elem.p2_id
        << ": " << v_map.at(std::to_string(elem.p2_id)); */
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
    auto id = std::stoi(inv_map.at(vertex));
    // BOOST_LOG_TRIVIAL(info) << "Vertex: " << vertex << " id: " << id;
    part_cell_map.insert(std::make_pair(vertex, part_cell_map.at(id)));
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

double GraphSQLite::get_mean()
{
  return boost::accumulators::mean(acc_total_domain);
}

Betweenness_Result GraphSQLite::get_result()
{
  std::map<std::string, std::map<std::string, double>> agg_per_cell{};
  for (auto &elem : decomp) {
    auto it = cell_agg.equal_range(elem.cellstr);
    if (it.first == it.second)
      continue;
    Acc cell_acc;
    for (auto i = it.first; i != it.second; ++i) {
      cell_acc((*i).second);
      // BOOST_LOG_TRIVIAL(info) << (*i).second;
    }
    agg_per_cell[elem.cellstr]["mean"] = accumulators::mean(cell_acc);
    agg_per_cell[elem.cellstr]["kur"] = accumulators::kurtosis(cell_acc);
    agg_per_cell[elem.cellstr]["var"] = accumulators::variance(cell_acc);
    agg_per_cell[elem.cellstr]["skew"] = accumulators::skewness(cell_acc);
  }

  Betweenness_Result res;
  res.aggr_per_cell = agg_per_cell;
  res.keys = keys;
  res.vals = vals;
  res.mean = accumulators::mean(acc_total_domain);
  res.ts = timestep;
  res.q_090 = get_percentile_vector(v_betweeness, 0.90);
  res.q_099 = get_percentile_vector(v_betweeness, 0.99);
  res.kur = accumulators::kurtosis(acc_total_domain);
  res.var = accumulators::variance(acc_total_domain);
  res.skew = accumulators::skewness(acc_total_domain);
  return res;
}
GraphSQLite::~GraphSQLite() {}