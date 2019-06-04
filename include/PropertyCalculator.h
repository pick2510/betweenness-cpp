//
// Created by std on 09.05.19.
//

#ifndef DEM_UTILS_PROPERTYCALCULATOR_H
#define DEM_UTILS_PROPERTYCALCULATOR_H

#include "grid.h"
#include <map>
#include <vector>

#include <Eigen/Dense>

enum m_contact_cols {
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

enum db_contact_cols {
  db_p1_id,
  db_p2_id,
  db_contact_overlap,
  db_cn_force_x,
  db_cn_force_y,
  db_cn_force_z,
  db_ct_force_x,
  db_ct_force_y,
  db_ct_force_z,
  db_cell_str,
  db_cell_ts,
  db_c_force_x,
  db_c_force_y,
  db_c_force_z,
};

constexpr int m_contact_cols_items = 18;

using d_mat = Eigen::MatrixXd;
using c_storage_index_t = decltype(indexContactStorage(""));
using p_storage_index_t = decltype(ParticleIndexStorage(""));
using tuple_contact_storage_t =
    decltype(std::declval<c_storage_index_t>().select(sqlite_orm::columns(
        &ContactColumns::p1_id, &ContactColumns::p2_id,
        &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
        &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
        &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
        &ContactColumns::ct_force_z, &ContactColumns::cellstr,
        &ContactColumns::ts, &ContactColumns::c_force_x,
        &ContactColumns::c_force_y, &ContactColumns::c_force_z)));

using tuple_particle_storage_t =
    decltype(std::declval<p_storage_index_t>().select(sqlite_orm::columns(
        &ParticleColumns::p_id, &ParticleColumns::p_vx, &ParticleColumns::p_vy,
        &ParticleColumns::p_vz, &ParticleColumns::p_coord,
        &ParticleColumns::p_disp_x, &ParticleColumns::p_disp_y,
        &ParticleColumns::p_disp_z, &ParticleColumns::p_ke_rot,
        &ParticleColumns::p_ke_tra, &ParticleColumns::p_ke_tra,
        &ParticleColumns::ts, &ParticleColumns::cellstr)));

class PropertyCalculator {
private:
  d_mat m_data;
  long ts;
  const std::map<int, double> &radius;
  tuple_contact_storage_t &contact_data;
  std::multimap<std::string, int> cellstridx{};
  aggregate_map_t aggregate_map{};
  int tuple_size;
  int vec_size;
  const decom_vec_storage_t &decomp_str;
  inline void calc_pot_energy(int &counter);

public:
  PropertyCalculator(const std::map<int, double> &radius,
                     tuple_contact_storage_t &contact_data,
                     tuple_particle_storage_t &particle_data,
                     const decom_vec_storage_t &decomp_str);
  virtual ~PropertyCalculator();
  void aggregate_per_cell();
  const aggregate_map_t &getAggregateMap() const;
};

#endif // DEM_UTILS_PROPERTYCALCULATOR_H
