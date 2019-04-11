#include "utils.h"

#include <unistd.h>
#include "natural_sort.hpp"
#include <map>
#include <glob.h>
#include <string.h>
#include <deque>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>

#include <boost/log/trivial.hpp>

#include "main.h"

std::vector<std::string> glob(const std::string &pattern)
{
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    vector<string> filenames;
    if (return_value == 0)
    {
        for (size_t i = 0; i < glob_result.gl_pathc; ++i)
        {
            filenames.push_back(string(glob_result.gl_pathv[i]));
        }
    }
    else if (return_value == 1 || return_value == 2)
    {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

std::deque<std::string> glob_deq(const std::string &pattern)
{
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    deque<string> filenames;
    if (return_value == 0)
    {
        for (size_t i = 0; i < glob_result.gl_pathc; ++i)
        {
            filenames.push_back(string(glob_result.gl_pathv[i]));
        }
    }
    else if (return_value == 1 || return_value == 2)
    {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


bool cmp_ts(const Result &a, const Result &b)
{
    return SI::natural::compare(std::to_string(a.ts), std::to_string(b.ts));
}

bool cmp_radii(const std::vector<std::string> &a, const std::vector<std::string> &b)
{
    int i_a = std::stoi(a[0]);
    int i_b = std::stoi(b[0]);
    return i_a < i_b;
}

std::vector<std::vector<std::string>> read_file(const std::string &file)
{
    const int skip_lines = 9;
    std::vector<std::vector<std::string>> file_content;
    std::ifstream radius_file;
    std::string line;
    radius_file.open(file);
    for (int i = 0; i < skip_lines; i++)
    {
        std::getline(radius_file, line);
    }
    while (!radius_file.eof())
    {
        std::vector<std::string> splitted_line;
        std::getline(radius_file, line);
        split_string(line, splitted_line);
        file_content.push_back(splitted_line);
    }
    radius_file.close();
    return file_content;
}

std::string trim(const std::string &s)
{
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) {
        return std::isspace(c);
    });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) {
                      return std::isspace(c);
                  })
                      .base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

std::map<std::string, int> get_vertice_map(const std::vector<std::vector<std::string>> &radiusfile)
{
    std::map<std::string, int> vertice_map;
    int i = 0;
    for (auto &elem : radiusfile)
    {
        vertice_map[elem[0]] = i++;
    }

    return vertice_map;
}

std::vector<double> get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile)
{
    int lt_size = std::stoi(radiusfile.back()[0]) + 1;
    std::vector<double> lookup_table(lt_size, 0);
    for (auto &elem : radiusfile)
    {
        lookup_table[std::stoi(elem[0])] = std::stod(elem[5]);
    }
    return lookup_table;
}

Config getCL(int &argc, char **argv)
{
    Config runningConfig;
    int opt;
    while ((opt = getopt(argc, argv, "i:o:s:")) != -1)
    {
        switch (opt)
        {
        case 'o':
            runningConfig.OutputPath = optarg;
            break;
        case 'i':
            runningConfig.InputPath = optarg;
            break;
        case 's':
            char *p = trimwhitespace(optarg);
            strncpy(runningConfig.sep, p, 1);
            runningConfig.sep[1] = '\0';
            break;
        }
    }
    if (runningConfig.OutputPath.empty() || runningConfig.InputPath.empty())
    {
        throw std::invalid_argument("Please use Input and Output Path");
    }
    return runningConfig;
}

void goto_line(std::ifstream &file, unsigned long n)
{
    std::string trashline;
    file.clear();
    file.seekg(0, std::ios::beg);

    for (unsigned long i = 0; i < n; i++)
    {
        std::getline(file, trashline);
    }
}

void write_ts_header(std::ofstream &out, Config &conf)
{

    out << "ts" << conf.sep << "mean" << conf.sep << "var" << conf.sep
        << "std" << conf.sep << "skew" << conf.sep
        << "kurtosis" << conf.sep << "q090" << conf.sep
        << "q099"
        << "\n";
}

void write_cent_header(std::ofstream &out, Config &conf)
{
    out << "particleid" << conf.sep << "centrality\n";
}

char *trimwhitespace(char *str)
{
    char *end;

    while (isspace((unsigned char)*str))
        str++;

    if (*str == 0)
        return str;

    // Trim trailing space
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end))
        end--;

    // Write new null terminator character
    end[1] = '\0';

    return str;
}

