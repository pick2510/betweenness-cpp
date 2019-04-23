#ifndef GRID_H
#define GRID_H
#include "data.h"
#include "sqlite_orm.h"

inline auto initStorage(const std::string &path)
{
  using namespace sqlite_orm;
  return make_storage(
      path,
      make_table(
          "ParticleContact", make_column("p1_x", &ContactColumns::p1_x),
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

#endif
