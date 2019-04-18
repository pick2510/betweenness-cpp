#ifndef UTIL_H
#define UTIL_H

#include "data.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <functional>
#include <iomanip>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> glob(const std::string &pattern);
std::deque<std::string> glob_deq(const std::string &pattern);
std::vector<std::vector<std::string>> read_file(const std::string &file);
bool cmp_radii(const std::vector<std::string> &a,
               const std::vector<std::string> &b);
std::string trim(const std::string &s);
std::map<std::string, int>
get_vertice_map(const std::vector<std::vector<std::string>> &radiusfile);
std::vector<double>
get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile);
Config getCL(int &argc, char **argv);
void goto_line(std::ifstream &file, unsigned long n);
bool cmp_ts(const Result &a, const Result &b);
void write_ts_header(std::ofstream &out, Config &conf);
void write_cent_header(std::ofstream &out, Config &conf);
void output_ts(std::ofstream &ts_mean_file, Config &runningConf,
               std::vector<Result> &results,
               std::map<int, std::string> &inv_vertice_map);
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
std::map<K, V> constructMap(std::vector<K> &keys, std::vector<V> &vals)
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

#endif