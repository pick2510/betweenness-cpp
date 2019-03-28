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



template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}



#endif