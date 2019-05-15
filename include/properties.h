//
// Created by std on 10.05.19.
//

#ifndef DEM_UTILS_PROPERTIES_H
#define DEM_UTILS_PROPERTIES_H
#include "data.h"
#include "grid.h"
#include <map>
#include <vector>
void calc_system_sum_write_file(std::vector<aggr_result_t> &results,
                                std::ofstream &f, Config &runningConf);
void write_cellstr_res(decom_vec_storage_t decomp_str,
                       std::vector<aggr_result_t> &results, Config &runningConf,
                       boost::filesystem::path &cellstr_path);
void initialize_output_files(const Config &runningConf,
                             decom_vec_storage_t decomp_str,
                             const boost::filesystem::path &system_path,
                             const boost::filesystem::path &cellstr_path);

void write_results(Config &runningConf, decom_vec_storage_t decomp_str,
                   std::vector<aggr_result_t> &results,
                   const boost::filesystem::path &system_path,
                   boost::filesystem::path &cellstr_path);

#endif // DEM_UTILS_PROPERTIES_H
