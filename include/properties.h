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
void create_cellstr_map(decom_vec_storage_t decomp_str,
                        std::vector<aggr_result_t> &results,
                        Config &runningConf,
                        boost::filesystem::path &cellstr_path);

#endif // DEM_UTILS_PROPERTIES_H
