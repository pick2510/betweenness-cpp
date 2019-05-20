//
// Created by std on 17.05.19.
//

#ifndef DEM_UTILS_VELOCITY_H
#define DEM_UTILS_VELOCITY_H
#include "grid.h"
#include <data.h>
#include <deque>
#include <set>
#include <vector>

std::set<int> select_particles(std::vector<int> &partids,
                               const Config &runningConf);
void process_ts(p_storage_index_t &particles, const std::vector<long> &ts,
                const std::set<int> &selected_particles,
                std::vector<p_velocity_t> &results, int ts_len);
void write_headers(const Config &runningConf,
                   const std::set<int> &selected_particles,
                   const boost::filesystem::path &out_path);

void write_results(const Config &runningConf,
                   const std::set<int> &selected_particles,
                   std::vector<p_velocity_t> &results,
                   const boost::filesystem::path &out_path, int ts_len);

#endif // DEM_UTILS_VELOCITY_H
