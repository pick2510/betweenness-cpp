//
// Created by std on 10.05.19.
//

#ifndef DEM_UTILS_PROPERTIES_H
#define DEM_UTILS_PROPERTIES_H
#include "data.h"
#include "grid.h"
#include <boost/mpi.hpp>
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

void submitCompleteWorld(const boost::mpi::communicator &world, int world_size,
                         int t_len, std::deque<long> &ts, c_storage_index_t &db,
                         std::vector<boost::mpi::request> &reqs_world,
                         std::vector<aggr_result_t> &results, long &v_index,
                         long &counter, p_storage_index_t &p_db);

void submitPieces(const boost::mpi::communicator &world, Config &runningConf,
                  int t_len, std::deque<long> &ts,
                  decom_vec_storage_t decomp_str,
                  std::vector<aggr_result_t> &results,
                  const boost::filesystem::path &system_path,
                  boost::filesystem::path &cellstr_path,
                  c_storage_index_t &c_db, int &chunk_len, long &v_index,
                  std::vector<boost::mpi::request> &reqs_world, bool &stop,
                  int world_size, long &counter, p_storage_index_t &p_db);

#endif // DEM_UTILS_PROPERTIES_H
