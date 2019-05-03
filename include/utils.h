#ifndef UTIL_H
#define UTIL_H

#include "INIReader.h"
#include "data.h"
#include <Eigen/Eigen>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <cassert>
#include <cmath>
#include <deque>
#include <functional>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> glob(const std::string &pattern);
std::deque<std::string> glob_deq(const std::string &pattern);
std::vector<std::vector<std::string>> read_radius_file(const std::string &file);
bool cmp_radii(const std::vector<std::string> &a,
               const std::vector<std::string> &b);
std::string trim(const std::string &s);
std::map<std::string, int>
get_vertice_map(const std::map<int, double> &radius_map);
std::map<std::string, int>
get_vertice_map(const std::vector<std::vector<std::string>> &radiusfile);
void get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile,
                      std::vector<double> &lookup_table,
                      std::map<int, double> &radius_map);
std::vector<double>
get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile);
std::vector<double> get_lookup_table(const std::map<int, double> &radius_map);
std::string getConfigPath(int &argc, char **argv);
Config getBetweennessConfigObj(INIReader &reader, const char *type);
INIReader parseConfigFile(const std::string &path);
void goto_line(std::ifstream &file, unsigned long n);
bool cmp_ts(const Result &a, const Result &b);
bool cmp_avg_map(const std::pair<std::string, double> &a,
                 const std::pair<std::string, double> &b);
void write_ts_header(std::ofstream &out, const Config &conf);
void write_cent_header(std::ofstream &out, const Config &conf);
void output_centrality_ts(std::ofstream &ts_mean_file,
                          const Config &runningConf,
                          const std::vector<Result> &results,
                          const std::map<int, std::string> &inv_vertice_map);
void output_particle_complete_ts(
    const Config &runningConf, const Eigen::Map<Eigen::MatrixXd> &mat,
    const std::map<int, std::string> &inv_vertice_map,
    const std::vector<long> &ts);
void output_particle_ts(const Config &runningConf,
                        const Eigen::Map<Eigen::MatrixXd> &mat,
                        const std::map<int, std::string> &inv_vertice_map,
                        const std::vector<long> &ts);
std::vector<std::string> getCartesianProduct(std::vector<int> &x,
                                             std::vector<int> &y,
                                             std::vector<int> &z);
void shuffleParticles(std::vector<int> &shuffleVec);
void LogConfig(Config &conf);

inline void check_path(const boost::filesystem::path &path);
Config getGridConfigObj(INIReader &reader);

char *trimwhitespace(char *str);

template <class Container>
void split_string(const std::string &str, Container &cont, char delim = ' ')
{
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delim)) {
    cont.push_back(token);
  }
}

template <typename K, typename V>
std::map<V, K> inverse_map(std::map<K, V> &map)
{
  std::map<V, K> inv;
  for (auto &kv : map) {
    inv[kv.second] = kv.first;
  }
  return inv;
}

template <typename K, typename V> std::vector<K> getKeys(std::map<K, V> &map)
{
  std::vector<K> keys;
  for (auto &kv : map) {
    keys.push_back(kv.first);
  }
  return keys;
}

template <typename K, typename V> std::vector<V> getVals(std::map<K, V> &map)
{
  std::vector<V> vals;
  for (auto &kv : map) {
    vals.push_back(kv.second);
  }
  return vals;
}

template <typename K, typename V>
std::map<K, V> constructMap(const std::vector<K> &keys,
                            const std::vector<V> &vals)
{
  std::map<K, V> v_map;
  std::transform(keys.begin(), keys.end(), vals.begin(),
                 std::inserter(v_map, v_map.end()),
                 [](K a, V b) { return std::make_pair(a, b); });
  return v_map;
}

template <typename K, typename V>
std::map<V, K> buildMap(std::vector<K> &keys, std::vector<V> &vals)
{
  assert(keys.size() == vals.size());
  std::map<V, K> kvmap;
  for (int i = 0; i < keys.size(); i++) {
    kvmap[keys[i]] = vals[i];
  }
  return kvmap;
}

template <typename T> std::string convFillString(T &num, int len)
{
  std::ostringstream oss;
  oss << std::setfill('0') << std::setw(len) << num;
  return oss.str();
}

template <typename K, typename V>
std::vector<std::pair<K, V>> get_vec_of_pairs(std::map<K, V> &map)
{
  std::vector<std::pair<K, V>> pair_vec;
  for (auto &kv : map) {
    pair_vec.push_back(std::make_pair(kv.first, kv.second));
  }
  return pair_vec;
}

template <typename T>
typename std::vector<T>::iterator
get_percentile_iterator(std::vector<T> &valVec, double percentile)
{
  auto n = std::ceil(percentile * valVec.size()) - 1;
  return valVec.begin() + n;
}

template <typename K, typename V>
bool cmp_second(const std::pair<K, V> &a, const std::pair<K, V> &b)
{
  return a.second < b.second;
}

#endif