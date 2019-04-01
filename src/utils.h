#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "data.h"

std::vector<std::string> glob(const std::string& pattern);
std::vector<std::vector<std::string>> read_file(const std::string& file);
bool cmp_radii(const std::vector<std::string> &a, const std::vector<std::string> &b);
std::string trim(const std::string &s);
std::map<std::string,int>get_vertice_map(const std::vector<std::vector<std::string>> &radiusfile);
std::vector<double> get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile);
Config getCL(int &argc, char **argv);
void goto_line(std::ifstream &file, unsigned long n);
bool cmp_ts(const Result &a, const Result &b);
void write_ts_header(std::ofstream &out);
void write_cent_header(std::ofstream &out);



template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}


template<typename K, typename V>
std::map<V,K> inverse_map(std::map<K,V> &map)
{
	std::map<V,K> inv;
	for (auto &kv: map){
        inv[kv.second] = kv.first;
    }
    return inv;
}



#endif