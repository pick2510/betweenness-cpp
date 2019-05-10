//
// Created by std on 09.05.19.
//

#ifndef DEM_UTILS_POTENTIAL_ENERGY_H
#define DEM_UTILS_POTENTIAL_ENERGY_H

#include <map>
#include <vector>

#include <Eigen/Dense>

enum m_cols {
  p1_r,
  p2_r,
  radstar,
  kn,
  kt,
  fnor,
  ftan,
  penor,
  petan,
  m_contact_overlap,
  m_cn_force_x,
  m_cn_force_y,
  m_cn_force_z,
  m_ct_force_x,
  m_ct_force_y,
  m_ct_force_z,
};

using d_mat = Eigen::MatrixXd;

class potential_energy {
  std::multimap<std::string, int> cellstridx;
  d_mat data;

public:
  potential_energy();
  virtual ~potential_energy();
};

#endif // DEM_UTILS_POTENTIAL_ENERGY_H
