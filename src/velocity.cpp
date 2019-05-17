//
// Created by std on 17.05.19.
//

#include "velocity.h"
#include "INIReader.h"
#include "PotentialEnergy.h"
#include "accumulator.h"
#include "data.h"
#include "grid.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <map>
#include <serialize_tuple.h>
#include <tuple>

int main(int argc, char *argv[])
{
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
  auto particles = ParticleIndexStorage(part_path);
  particles.sync_schema();
  auto radius_storage = initRadstorage(path);
  for (auto &elem : radius_storage.iterate<radius>()) {
    radius_map[elem.id] = elem.rad;
    partids.push_back(elem.id);
  }
  auto ts_list = stor.get_all<ts_column>();
  std::deque<long> ts;
  for (const auto &elem : ts_list) {
    ts.push_back(elem.ts);
  }
  shuffleParticles(partids);
  assert(runningConf.randomly_selected <= partids.size());
  for (auto i = partids.begin();
       i < partids.begin() + runningConf.randomly_selected; ++i) {
    BOOST_LOG_TRIVIAL(info) << *i;
  }
}