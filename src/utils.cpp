#include "utils.h"

#include "natural_sort.hpp"
#include <cstring>
#include <deque>
#include <fstream>
#include <glob.h>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string.h>
#include <unistd.h>

#include "INIReader.h"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>

namespace fs = boost::filesystem;

std::vector<std::string> glob(const std::string &pattern)
{
  using namespace std;

  // glob struct resides on the stack
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));

  // do the glob operation
  int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  vector<string> filenames;
  if (return_value == 0) {
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
      filenames.push_back(string(glob_result.gl_pathv[i]));
    }
  }
  else if (return_value == 1 || return_value == 2) {
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
  if (return_value == 0) {
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
      filenames.push_back(string(glob_result.gl_pathv[i]));
    }
  }
  else if (return_value == 1 || return_value == 2) {
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

bool cmp_radii(const std::vector<std::string> &a,
               const std::vector<std::string> &b)
{
  int i_a = std::stoi(a[0]);
  int i_b = std::stoi(b[0]);
  return i_a < i_b;
}

std::vector<std::vector<std::string>> read_radius_file(const std::string &file)
{
  const int skip_lines = 9;
  std::vector<std::vector<std::string>> file_content;
  std::ifstream radius_file;
  std::string line;
  radius_file.open(file);
  for (int i = 0; i < skip_lines; i++) {
    std::getline(radius_file, line);
  }
  while (!radius_file.eof()) {
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
  auto wsfront = std::find_if_not(s.begin(), s.end(),
                                  [](int c) { return std::isspace(c); });
  auto wsback = std::find_if_not(s.rbegin(), s.rend(),
                                 [](int c) { return std::isspace(c); })
                    .base();
  return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

std::map<std::string, int>
get_vertice_map(const std::vector<std::vector<std::string>> &radiusfile)
{
  std::map<std::string, int> vertice_map;
  int i = 0;
  for (auto &elem : radiusfile) {
    vertice_map[elem[0]] = i++;
  }

  return vertice_map;
}

std::map<std::string, int>
get_vertice_map(const std::map<int, double> &radius_map)
{
  std::map<std::string, int> vertice_map;
  int i = 0;
  for (auto &elem : radius_map) {
    vertice_map[std::to_string(elem.first)] = i++;
  }

  return vertice_map;
}

std::vector<double> get_lookup_table(const std::map<int, double> &radius_map)
{
  std::vector<double> lookup_table(radius_map.rbegin()->first + 1, 0);
  for (const auto &elem : radius_map) {
    lookup_table[elem.first] = elem.second;
  }
  return lookup_table;
}

void get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile,
                      std::vector<double> &lookup_table,
                      std::map<int, double> &radius_map)
{
  auto lt_size = std::stoi(radiusfile.back()[0]) + 1;
  lookup_table.resize(lt_size);
  for (auto &elem : radiusfile) {
    lookup_table[std::stoi(elem[0])] = std::stod(elem[5]);
    radius_map.insert(
        std::make_pair<int, double>(std::stoi(elem[0]), std::stod(elem[5])));
  }
}

std::vector<double>
get_lookup_table(const std::vector<std::vector<std::string>> &radiusfile)
{
  std::vector<double> lookup_table;
  auto lt_size = std::stoi(radiusfile.back()[0]) + 1;
  lookup_table.resize(lt_size);
  for (auto &elem : radiusfile) {
    lookup_table[std::stoi(elem[0])] = std::stod(elem[5]);
  }
  return lookup_table;
}

std::string getConfigPath(int &argc, char **argv)
{
  std::string ConfigPath{};
  int opt;
  while ((opt = getopt(argc, argv, "c:")) != -1) {
    switch (opt) {
    case 'c':
      ConfigPath = optarg;
      break;
    }
  }
  if (ConfigPath.empty()) {
    throw std::invalid_argument("Please -c Path to configfile");
  }
  return ConfigPath;
}

INIReader parseConfigFile(const std::string &path)
{
  fs::path p(path);
  if (!fs::exists(p)) {
    throw std::invalid_argument("Configfile doesn't exist");
  }
  INIReader reader(p.string());
  if (reader.ParseError() != 0) {
    throw std::invalid_argument("Configfile " + p.string() +
                                " couldn't be loaded");
  }
  return reader;
}

void goto_line(std::ifstream &file, unsigned long n)
{
  std::string trashline;
  file.clear();
  file.seekg(0, std::ios::beg);

  for (unsigned long i = 0; i < n; i++) {
    std::getline(file, trashline);
  }
}

void write_ts_header(std::ofstream &out, const Config &conf)
{

  out << "ts" << conf.sep << "mean" << conf.sep << "var" << conf.sep << "std"
      << conf.sep << "skew" << conf.sep << "kurtosis" << conf.sep << "q090"
      << conf.sep << "q099"
      << "\n";
}

void write_cent_header(std::ofstream &out, const Config &conf)
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

void check_path(const boost::filesystem::path &path)
{
  if (!fs::exists(path)) {
    fs::create_directories(path);
  }
}

void output_centrality_ts(std::ofstream &ts_mean_file,
                          const Config &runningConf,
                          const std::vector<Result> &results,
                          const std::map<int, std::string> &inv_vertice_map)
{
  auto p = fs::path(runningConf.OutputPath + "/" + ts_centrality_path);
  check_path(p);
  for (auto &v : results) {
    ts_mean_file << convFillString(v.ts, 9) << runningConf.sep
                 << std::to_string(v.mean) << runningConf.sep
                 << std::to_string(v.var) << runningConf.sep
                 << std::to_string(std::sqrt(v.var)) << runningConf.sep
                 << std::to_string(v.skew) << runningConf.sep
                 << std::to_string(v.kur) << runningConf.sep
                 << std::to_string(v.q_090) << runningConf.sep
                 << std::to_string(v.q_099) << "\n";
    std::ofstream ts_file(p.string() + "/centrality_" +
                          convFillString(v.ts, 9) + ".csv");
    write_cent_header(ts_file, runningConf);
    ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto b_centrality = constructMap(v.keys, v.vals);
    for (auto &kv : b_centrality) {
      ts_file << inv_vertice_map.at(kv.first) << runningConf.sep << kv.second
              << "\n";
    }
    ts_file.close();
  }
}

void output_particle_complete_ts(
    const Config &runningConf, const Eigen::Map<Eigen::MatrixXd> &mat,
    const std::map<int, std::string> &inv_vertice_map,
    const std::vector<long> &ts)
{
  auto p = fs::path(runningConf.OutputPath + "/" + ts_particle_path);
  check_path(p);
  auto le = mat.cols();
  auto ts_len = ts.size();
  for (unsigned int i = 0; i < le; i++) {
    auto col = mat.col(i);
    std::ofstream ts_file(p.string() + "/" +
                          convFillString(inv_vertice_map.at(i), 5) + ".csv");
    ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    ts_file << "ts" << runningConf.sep << "centrality\n";
    for (unsigned int j = 0; j < ts_len; j++) {
      ts_file << ts[j] << runningConf.sep << col(j) << "\n";
    }
    ts_file.close();
  }
}

void output_particle_ts(const Config &runningConf,
                        const Eigen::Map<Eigen::MatrixXd> &mat,
                        const std::map<int, std::string> &inv_vertice_map,
                        const std::vector<long> &ts)
{
  auto p = fs::path(runningConf.OutputPath + "/" + ts_particle_path);
  auto p_q = p / "quantiles";
  auto p_r = p / "random";
  check_path(p_q);
  check_path(p_r);
  auto le = mat.cols();
  auto ts_len = ts.size();
  std::vector<std::pair<unsigned int, double>> averages_pair;
  std::vector<int> ids;
  auto averages_mat = mat.colwise().sum() / mat.rows();
  for (unsigned int i = 0; i < le; i++) {
    auto col = averages_mat.col(i);
    averages_pair.push_back(std::make_pair(i, col.value()));
    ids.push_back(i);
  }
  std::sort(averages_pair.begin(), averages_pair.end(),
            cmp_second<unsigned int, double>);
  auto n =
      get_percentile_iterator(averages_pair, runningConf.output_percentile);
  for (auto i = n; i != averages_pair.end(); ++i) {
    auto col = mat.col((*i).first);
    std::ofstream ts_file(p_q.string() + "/" +
                          convFillString(inv_vertice_map.at((*i).first), 5) +
                          ".csv");
    ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    ts_file << "ts" << runningConf.sep << "centrality\n";
    for (unsigned int j = 0; j < ts_len; j++) {
      ts_file << ts[j] << runningConf.sep << col(j) << "\n";
    }
    ts_file.close();
  }
  shuffleParticles(ids);
  assert(runningConf.randomly_selected <= ids.size());

  for (auto i = ids.begin(); i < ids.begin() + runningConf.randomly_selected;
       ++i) {
    auto col = mat.col(*i);
    std::ofstream ts_file(p_r.string() + "/" +
                          convFillString(inv_vertice_map.at(*i), 5) + ".csv");
    ts_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    ts_file << "ts" << runningConf.sep << "centrality\n";
    for (unsigned int j = 0; j < ts_len; j++) {
      ts_file << ts[j] << runningConf.sep << col(j) << "\n";
    }

    ts_file.close();
  }
}

void shuffleParticles(std::vector<int> &shuffleVec)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  auto currentIndexCounter = shuffleVec.size();
  for (auto iter = shuffleVec.rbegin(); iter != shuffleVec.rend();
       ++iter, --currentIndexCounter) {
    std::uniform_int_distribution<> dis(0, currentIndexCounter);
    const int randomIndex = dis(mt);

    if (*iter != shuffleVec.at(randomIndex)) {
      std::swap(shuffleVec.at(randomIndex), *iter);
    }
  }
}

Config getGridConfigObj(INIReader &reader)
{
  Config conf;
  conf.InputPath = reader.Get("grid", "inputPath", "");
  conf.OutputPath = reader.Get("grid", "outputPath", "");
  conf.x_cells = static_cast<int>(reader.GetInteger("grid", "x_cells", 0));
  conf.y_cells = static_cast<int>(reader.GetInteger("grid", "y_cells", 0));
  conf.z_cells = static_cast<int>(reader.GetInteger("grid", "z_cells", 0));
  conf.domainsize_x = reader.GetReal("grid", "domainsize_x", 0.0);
  conf.domainsize_y = reader.GetReal("grid", "domainsize_y", 0.0);
  conf.domainsize_z = reader.GetReal("grid", "domainsize_z", 0.0);
  conf.contact_filename = reader.Get("grid", "ContactfileName", "");
  conf.particle_filename = reader.Get("grid", "ParticlefileName", "");
  return conf;
}

Config getBetweennessConfigObj(INIReader &reader, const char *type)
{
  Config conf;
  conf.InputPath = reader.Get(type, "inputPath", "");
  conf.OutputPath = reader.Get(type, "outputPath", "");
  conf.output_percentile = reader.GetReal(type, "outputPercentile", 0.9);
  conf.randomly_selected =
      static_cast<int>(reader.GetInteger(type, "randomlySelected", 20));
  conf.contact_filename = reader.Get(type, "ContactfileName", "");
  return conf;
}

Config getPropertiesConfigObj(INIReader &reader){
    Config conf;
    conf.InputPath = reader.Get("properties", "inputPath", "");
    conf.OutputPath = reader.Get("properties", "outputPath", "");
    conf.contact_filename = reader.Get("properties", "ContactfileName", "");
    conf.particle_filename = reader.Get("properties", "ParticlefileName", "");
    conf.spinup_time = reader.GetInteger("properties", "spinupTime", 0);
    return conf;
}

std::vector<std::string> getCartesianProduct(std::vector<int> &x,
                                             std::vector<int> &y,
                                             std::vector<int> &z)
{
  std::vector<std::string> cellset;
  for (auto &elem_x : x) {
    for (auto &elem_y : y) {
      for (auto &elem_z : z) {
        std::stringstream ss;
        ss << "cell_" << std::setfill('0') << std::setw(2) << elem_x << "_"
           << std::setfill('0') << std::setw(2) << elem_y << "_"
           << std::setfill('0') << std::setw(2) << elem_z;
        cellset.push_back(ss.str());
      }
    }
  }
  return cellset;
}

void LogConfig(Config &conf)
{
  BOOST_LOG_TRIVIAL(info) << "Using Input: " << conf.InputPath;
  BOOST_LOG_TRIVIAL(info) << "Using Output: " << conf.OutputPath;
  BOOST_LOG_TRIVIAL(info) << "x_cells: " << conf.x_cells;
  BOOST_LOG_TRIVIAL(info) << "y_cells: " << conf.y_cells;
  BOOST_LOG_TRIVIAL(info) << "z_cells: " << conf.z_cells;
  BOOST_LOG_TRIVIAL(info) << "Domain size: " << conf.domainsize_x << "x"
                          << conf.domainsize_y << "x" << conf.domainsize_z;
}
