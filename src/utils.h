#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <functional>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "main.h"

std::vector<std::string> glob(const std::string& pattern);
std::vector<std::vector<std::string>> read_file(const std::string& file);
bool cmp_radii(const std::vector<std::string> &a, const std::vector<std::string> &b);
std::string trim(const std::string &s);
std::vector<double> get_lookup_table(std::vector<std::vector<std::string>> &radiusfile);
Config getCL(int &argc, char **argv);



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