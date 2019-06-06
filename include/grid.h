#ifndef GRID_H
#define GRID_H
#include "data.h"
#include "decomposition.h"
#include "dumpfile.h"
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "sqlite_orm.h"
#include <atomic>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include <vector>

void LogConfig(Config &conf);

void processContacts(std::deque<std::string> &chain_file_list,
                     std::map<int, double> &radius_map, Decomposition &decomp,
                     std::vector<long> &ts_vec, std::size_t chain_size,
                     std::atomic<long> &index,
                     std::vector<std::vector<ContactColumns>> &chunk_res);
void processParticles(std::vector<std::string> &radius_file_list,
                      std::map<int, double> &radius_map, Decomposition &decomp,
                      std::size_t part_size, std::atomic<long> &index,
                      std::vector<std::vector<ParticleColumns>> &chunk_res);

inline auto initContactStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path,
      make_table(
          "ParticleContact",
          make_column("id", &ContactColumns::id, primary_key()),
          make_column("p1_x", &ContactColumns::p1_x),
          make_column("p1_y", &ContactColumns::p1_y),
          make_column("p1_z", &ContactColumns::p1_z),
          make_column("p2_x", &ContactColumns::p2_x),
          make_column("p2_y", &ContactColumns::p2_y),
          make_column("p2_z", &ContactColumns::p2_z),
          make_column("p1_vx", &ContactColumns::p1_vx),
          make_column("p1_vy", &ContactColumns::p1_vy),
          make_column("p1_vz", &ContactColumns::p1_vz),
          make_column("p2_vx", &ContactColumns::p2_vx),
          make_column("p2_vy", &ContactColumns::p2_vy),
          make_column("p2_vz", &ContactColumns::p2_vz),
          make_column("p1_id", &ContactColumns::p1_id),
          make_column("p2_id", &ContactColumns::p2_id),
          make_column("is_periodic", &ContactColumns::is_periodic),
          make_column("c_force_x", &ContactColumns::c_force_x),
          make_column("c_force_y", &ContactColumns::c_force_y),
          make_column("c_force_z", &ContactColumns::c_force_z),
          make_column("cn_force_x", &ContactColumns::cn_force_x),
          make_column("cn_force_y", &ContactColumns::cn_force_y),
          make_column("cn_force_z", &ContactColumns::cn_force_z),
          make_column("ct_force_x", &ContactColumns::ct_force_x),
          make_column("ct_force_y", &ContactColumns::ct_force_y),
          make_column("ct_force_z", &ContactColumns::ct_force_z),
          make_column("c_torque_x", &ContactColumns::c_torque_x),
          make_column("c_torque_y", &ContactColumns::c_torque_y),
          make_column("c_torque_z", &ContactColumns::c_torque_z),
          make_column("disp_x", &ContactColumns::disp_x),
          make_column("disp_y", &ContactColumns::disp_y),
          make_column("disp_z", &ContactColumns::disp_z),
          make_column("contact_area", &ContactColumns::contact_area),
          make_column("contact_overlap", &ContactColumns::contact_overlap),
          make_column("sliding_contact", &ContactColumns::sliding_contact),
          make_column("cell_x", &ContactColumns::cell_x),
          make_column("cell_y", &ContactColumns::cell_y),
          make_column("cell_z", &ContactColumns::cell_z),
          make_column("ts", &ContactColumns::ts),
          make_column("cellstr", &ContactColumns::cellstr)));
}

inline auto initTSStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(path,
                      make_table("TS", make_column("ts", &ts_column::ts)));
}

inline auto initRadstorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(path, make_table("Radius",
                                       make_column("particleid", &radius::id),
                                       make_column("radius", &radius::rad)));
}

inline auto initDecompstorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path, make_table("Decomposition",
                       make_column("cellstr", &decomp_table::cellstr)));
}

inline auto indexContactStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path, make_index("idx_ts_cellstr", &ContactColumns::ts),
      make_table(
          "ParticleContact",
          make_column("id", &ContactColumns::id, primary_key()),
          make_column("p1_x", &ContactColumns::p1_x),
          make_column("p1_y", &ContactColumns::p1_y),
          make_column("p1_z", &ContactColumns::p1_z),
          make_column("p2_x", &ContactColumns::p2_x),
          make_column("p2_y", &ContactColumns::p2_y),
          make_column("p2_z", &ContactColumns::p2_z),
          make_column("p1_vx", &ContactColumns::p1_vx),
          make_column("p1_vy", &ContactColumns::p1_vy),
          make_column("p1_vz", &ContactColumns::p1_vz),
          make_column("p2_vx", &ContactColumns::p2_vx),
          make_column("p2_vy", &ContactColumns::p2_vy),
          make_column("p2_vz", &ContactColumns::p2_vz),
          make_column("p1_id", &ContactColumns::p1_id),
          make_column("p2_id", &ContactColumns::p2_id),
          make_column("is_periodic", &ContactColumns::is_periodic),
          make_column("c_force_x", &ContactColumns::c_force_x),
          make_column("c_force_y", &ContactColumns::c_force_y),
          make_column("c_force_z", &ContactColumns::c_force_z),
          make_column("cn_force_x", &ContactColumns::cn_force_x),
          make_column("cn_force_y", &ContactColumns::cn_force_y),
          make_column("cn_force_z", &ContactColumns::cn_force_z),
          make_column("ct_force_x", &ContactColumns::ct_force_x),
          make_column("ct_force_y", &ContactColumns::ct_force_y),
          make_column("ct_force_z", &ContactColumns::ct_force_z),
          make_column("c_torque_x", &ContactColumns::c_torque_x),
          make_column("c_torque_y", &ContactColumns::c_torque_y),
          make_column("c_torque_z", &ContactColumns::c_torque_z),
          make_column("disp_x", &ContactColumns::disp_x),
          make_column("disp_y", &ContactColumns::disp_y),
          make_column("disp_z", &ContactColumns::disp_z),
          make_column("contact_area", &ContactColumns::contact_area),
          make_column("contact_overlap", &ContactColumns::contact_overlap),
          make_column("sliding_contact", &ContactColumns::sliding_contact),
          make_column("cell_x", &ContactColumns::cell_x),
          make_column("cell_y", &ContactColumns::cell_y),
          make_column("cell_z", &ContactColumns::cell_z),
          make_column("ts", &ContactColumns::ts),
          make_column("cellstr", &ContactColumns::cellstr)));
}

inline auto ParticleStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path, make_index("idx_ts_cellstr", &ParticleColumns::ts),
      make_table("Particles",
                 make_column("id", &ParticleColumns::id, primary_key()),
                 make_column("p_id", &ParticleColumns::p_id),
                 make_column("p_type", &ParticleColumns::p_type),
                 make_column("p_x", &ParticleColumns::p_x),
                 make_column("p_y", &ParticleColumns::p_y),
                 make_column("p_z", &ParticleColumns::p_z),
                 make_column("p_rad", &ParticleColumns::p_rad),
                 make_column("p_vx", &ParticleColumns::p_vx),
                 make_column("p_vy", &ParticleColumns::p_vy),
                 make_column("p_vz", &ParticleColumns::p_vz),
                 make_column("p_fx", &ParticleColumns::p_fx),
                 make_column("p_fy", &ParticleColumns::p_fy),
                 make_column("p_fz", &ParticleColumns::p_fz),
                 make_column("p_omegax", &ParticleColumns::p_omegax),
                 make_column("p_omegay", &ParticleColumns::p_omegay),
                 make_column("p_omegaz", &ParticleColumns::p_omegaz),
                 make_column("p_coord", &ParticleColumns::p_coord),
                 make_column("p_disp_x", &ParticleColumns::p_disp_x),
                 make_column("p_disp_y", &ParticleColumns::p_disp_y),
                 make_column("p_disp_z", &ParticleColumns::p_disp_z),
                 make_column("p_disp_mag", &ParticleColumns::p_disp_mag),
                 make_column("p_ke_rot", &ParticleColumns::p_ke_rot),
                 make_column("p_ke_tra", &ParticleColumns::p_ke_tra),
                 make_column("cell_x", &ParticleColumns::cell_x),
                 make_column("cell_y", &ParticleColumns::cell_y),
                 make_column("cell_z", &ParticleColumns::cell_z),
                 make_column("ts", &ParticleColumns::ts),
                 make_column("cellstr", &ParticleColumns::cellstr)));
}

inline auto ParticleIndexStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path, make_table("Particles",
                       make_column("id", &ParticleColumns::id, primary_key()),
                       make_column("p_id", &ParticleColumns::p_id),
                       make_column("p_type", &ParticleColumns::p_type),
                       make_column("p_x", &ParticleColumns::p_x),
                       make_column("p_y", &ParticleColumns::p_y),
                       make_column("p_z", &ParticleColumns::p_z),
                       make_column("p_rad", &ParticleColumns::p_rad),
                       make_column("p_vx", &ParticleColumns::p_vx),
                       make_column("p_vy", &ParticleColumns::p_vy),
                       make_column("p_vz", &ParticleColumns::p_vz),
                       make_column("p_fx", &ParticleColumns::p_fx),
                       make_column("p_fy", &ParticleColumns::p_fy),
                       make_column("p_fz", &ParticleColumns::p_fz),
                       make_column("p_omegax", &ParticleColumns::p_omegax),
                       make_column("p_omegay", &ParticleColumns::p_omegay),
                       make_column("p_omegaz", &ParticleColumns::p_omegaz),
                       make_column("p_coord", &ParticleColumns::p_coord),
                       make_column("p_disp_x", &ParticleColumns::p_disp_x),
                       make_column("p_disp_y", &ParticleColumns::p_disp_y),
                       make_column("p_disp_z", &ParticleColumns::p_disp_z),
                       make_column("p_disp_mag", &ParticleColumns::p_disp_mag),
                       make_column("p_ke_rot", &ParticleColumns::p_ke_rot),
                       make_column("p_ke_tra", &ParticleColumns::p_ke_tra),
                       make_column("cell_x", &ParticleColumns::cell_x),
                       make_column("cell_y", &ParticleColumns::cell_y),
                       make_column("cell_z", &ParticleColumns::cell_z),
                       make_column("ts", &ParticleColumns::ts),
                       make_column("cellstr", &ParticleColumns::cellstr)));
}

using c_storage_t = decltype(initContactStorage(""));
using p_storage_t = decltype(ParticleStorage(""));
using c_storage_index_t = decltype(indexContactStorage(""));
using p_storage_index_t = decltype(ParticleIndexStorage(""));
using tuple_contact_storage_t =
    decltype(std::declval<c_storage_index_t>().select(sqlite_orm::columns(
        &ContactColumns::p1_id, &ContactColumns::p2_id,
        &ContactColumns::contact_overlap, &ContactColumns::cn_force_x,
        &ContactColumns::cn_force_y, &ContactColumns::cn_force_z,
        &ContactColumns::ct_force_x, &ContactColumns::ct_force_y,
        &ContactColumns::ct_force_z, &ContactColumns::c_force_x,
        &ContactColumns::c_force_y, &ContactColumns::c_force_z,
        &ContactColumns::sliding_contact, &ContactColumns::cellstr,
        &ContactColumns::ts)));

using tuple_particle_storage_t =
    decltype(std::declval<p_storage_index_t>().select(sqlite_orm::columns(
        &ParticleColumns::p_id, &ParticleColumns::p_vx, &ParticleColumns::p_vy,
        &ParticleColumns::p_vz, &ParticleColumns::p_coord,
        &ParticleColumns::p_disp_x, &ParticleColumns::p_disp_y,
        &ParticleColumns::p_disp_z, &ParticleColumns::p_ke_rot,
        &ParticleColumns::p_ke_tra, &ParticleColumns::p_omegax,
        &ParticleColumns::p_omegay, &ParticleColumns::p_omegaz,
        &ParticleColumns::cellstr, &ParticleColumns::ts)));

using decom_storage_t = decltype(initDecompstorage(""));
using decom_vec_storage_t =
    decltype(std::declval<decom_storage_t>().get_all<decomp_table>());

#endif
