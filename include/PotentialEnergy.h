//
// Created by std on 09.05.19.
//

#ifndef DEM_UTILS_POTENTIALENERGY_H
#define DEM_UTILS_POTENTIALENERGY_H

#include "grid.h"
#include <map>
#include <vector>

#include <Eigen/Dense>

enum m_cols {
  m_p1_id,
  m_p2_id,
  m_contact_overlap,
  m_cn_force_x,
  m_cn_force_y,
  m_cn_force_z,
  m_ct_force_x,
  m_ct_force_y,
  m_ct_force_z,
  m_p1_r,
  m_p2_r,
  m_radstar,
  m_kn,
  m_kt,
  m_fnor,
  m_ftan,
  m_penor,
  m_petan,
};

constexpr int m_cols_items = 18;

using d_mat = Eigen::MatrixXd;
using c_storage_index_t = decltype(indexContactStorage(""));
using tuple_storage_t =
    decltype(std::declval<c_storage_index_t>().select(sqlite_orm::columns(
        &ContactColumns::p1_id, &ContactColumns::p2_id,
        &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
        &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
        &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
        &ContactColumns::ct_force_z, &ContactColumns::cellstr,
        &ContactColumns::ts)));

class PotentialEnergy {
private:
  d_mat m_data;
  long ts;
  const std::map<int, double> &radius;
  tuple_storage_t &data;
  std::multimap<std::string, int> cellstridx{};
  std::map<std::string, double> aggregate{};
  int tuple_size;
  int vec_size;
  const decom_vec_storage_t &decomp_str;
  inline void calc_pot_energy(int &counter);

public:
  PotentialEnergy(const std::map<int, double> &radius, tuple_storage_t &data,
                  const decom_vec_storage_t &decomp_str);
  virtual ~PotentialEnergy();
  void aggregate_per_cell();
};

#endif // DEM_UTILS_POTENTIALENERGY_H
