//
// Created by std on 17.05.19.
//

#include "velocity.h"
#include "INIReader.h"
#include "PropertyCalculator.h"
#include "accumulator.h"
#include "data.h"
#include "grid.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iterator>
#include <map>
#include <serialize_tuple.h>
#include <tuple>

int main(int argc, char *argv[])
{
  using namespace sqlite_orm;
  std::map<int, double> radius_map;
  char hostname[HOSTNAME_LEN]{};
  std::string configPath{};
  std::vector<int> partids;
  INIReader reader;
  gethostname(hostname, HOSTNAME_LEN);
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Velocity Extraction of DEM Data";
  BOOST_LOG_TRIVIAL(info) << "by Dominik Strebel (2019)";
  BOOST_LOG_TRIVIAL(info) << "****************************************";
  BOOST_LOG_TRIVIAL(info) << "Hostname MASTER: " << hostname;
  try {
    configPath = getConfigPath(argc, argv);
  }
  catch (std::invalid_argument &e) {
    BOOST_LOG_TRIVIAL(error) << e.what();
    exit(EXIT_FAILURE);
  }
  try {
    reader = parseConfigFile(configPath);
  }
  catch (std::invalid_argument &e) {
    BOOST_LOG_TRIVIAL(error) << e.what();
    exit(EXIT_FAILURE);
  }
  auto runningConf = getVelocityConfigObj(reader);

  auto path =
      std::string{runningConf.InputPath + "/" + runningConf.contact_filename};
  auto part_path =
      std::string{runningConf.InputPath + "/" + runningConf.particle_filename};
  auto stor = initTSStorage(path);
  auto particle_storage = ParticleIndexStorage(part_path);
  auto radius_storage = initRadstorage(path);
  for (auto &elem : radius_storage.iterate<radius>()) {
    radius_map[elem.id] = elem.rad;
    partids.push_back(elem.id);
  }
  auto ts_list = stor.get_all<ts_column>();
  std::vector<long> ts;
  for (const auto &elem : ts_list) {
    ts.push_back(elem.ts);
  }
  auto ts_len = ts.size();
  assert(runningConf.randomly_selected <= partids.size());
  auto selected_particles = select_particles(partids, runningConf);
  std::vector<p_velocity_t> results{};
  process_ts(particle_storage, ts, selected_particles, results, ts_len);
  auto out_path = boost::filesystem::path(runningConf.OutputPath +
                                          std::string{"/velocities"});
  BOOST_LOG_TRIVIAL(info) << out_path.string();
  check_path(out_path);
  write_headers(runningConf, selected_particles, out_path);
  write_results(runningConf, selected_particles, results, out_path, 0);
}
void write_results(const Config &runningConf,
                   const std::set<int> &selected_particles,
                   std::vector<p_velocity_t> &results,
                   const boost::filesystem::path &out_path, int ts_len)
{
  std::map<int, std::ofstream> f_streams;
  for (const auto &part : selected_particles) {
    f_streams[part] = std::ofstream(
        out_path.string() + "/" + std::to_string(part) + ".csv", std::ios::app);
  }
  int counter = 0;
  for (auto &res : results) {

    for (const auto &part : selected_particles) {
      f_streams.at(part) << res.ts << runningConf.sep << res.vel[part]["p_vx"]
                         << runningConf.sep << res.vel[part]["p_vy"]
                         << runningConf.sep << res.vel[part]["p_vz"]
                         << runningConf.sep << res.vel[part]["mag_v"] << "\n";
    }
    BOOST_LOG_TRIVIAL(info)
        << static_cast<double>(counter) / ts_len * 100.0 << "% ";
    counter++;
  }
}
void write_headers(const Config &runningConf,
                   const std::set<int> &selected_particles,
                   const boost::filesystem::path &out_path)
{
  for (const auto &elem : selected_particles) {
    std::ofstream f(out_path.string() + "/" + std::to_string(elem) + ".csv");
    f << "ts" << runningConf.sep << "vx" << runningConf.sep << "vy"
      << runningConf.sep << "vz" << runningConf.sep << "mag(v)\n";
  }
}
void process_ts(p_storage_index_t &particles, const std::vector<long> &ts,
                const std::set<int> &selected_particles,
                std::vector<p_velocity_t> &results, int ts_len)
{
  using namespace sqlite_orm;
  int counter = 0;
  for (auto &elem : ts) {
    p_velocity_t p_velocity;
    p_velocity.ts = elem;
    std::map<int, std::map<std::string, double>> p_velocities_m{};
    auto timestep =
        particles.select(columns(&ParticleColumns::ts, &ParticleColumns::p_id,
                                 &ParticleColumns::p_vx, &ParticleColumns::p_vy,
                                 &ParticleColumns::p_vz),
                         where(c(&ParticleColumns::ts) == elem));
    BOOST_LOG_TRIVIAL(info) << "Timestep: " << p_velocity.ts;
    for (auto &particle : timestep) {
      auto it = selected_particles.find(std::get<1>(particle));
      if (it == selected_particles.end())
        continue;
      p_velocities_m[*it]["p_vx"] = std::get<2>(particle);
      p_velocities_m[*it]["p_vy"] = std::get<3>(particle);
      p_velocities_m[*it]["p_vz"] = std::get<4>(particle);
      p_velocities_m[*it]["mag_v"] =
          std::sqrt(std::pow(std::get<2>(particle), 2) +
                    std::pow(std::get<3>(particle), 2) +
                    std::pow(std::get<4>(particle), 2));
    }
    p_velocity.vel = p_velocities_m;
    results.push_back(p_velocity);
    BOOST_LOG_TRIVIAL(info)
        << static_cast<double>(counter) / ts_len * 100.0 << "% ";
    counter++;
  }
}
std::set<int> select_particles(std::vector<int> &partids,
                               const Config &runningConf)
{
  shuffleParticles(partids);
  std::set<int> selected_particles;
  std::copy(partids.begin(), partids.begin() + runningConf.randomly_selected,
            std::inserter(selected_particles, selected_particles.begin()));
  BOOST_LOG_TRIVIAL(info) << "Using particles: ";
  for (const auto &elem : selected_particles) {
    BOOST_LOG_TRIVIAL(info) << elem;
  }
  return selected_particles;
}