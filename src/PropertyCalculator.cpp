//
// Created by std on 09.05.19.
//

#include "PropertyCalculator.h"
#include "accumulator.h"
#include <Eigen/Dense>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <tuple>

PropertyCalculator::PropertyCalculator(const std::map<int, double> &radius,
                                       tuple_contact_storage_t &contact_data,
                                       tuple_particle_storage_t &particle_data,
                                       const decom_vec_storage_t &decomp_str)
    : radius(radius), contact_data(contact_data),
      vec_c_size(contact_data.size()), vec_p_size(particle_data.size()),
      decomp_str(decomp_str)
{
  // tuple_size =
  // std::tuple_size<std::decay_t<decltype(contact_data[0])>>::value;
  m_c_data = d_mat(vec_c_size, m_contact_cols_items);
  ts = std::get<14>(contact_data[0]);
  int counter = 0;
  for (const auto &elem : contact_data) {
    c_cellstridx.insert(
        std::pair<std::string, int>(std::get<13>(elem), counter));
    m_c_data(counter, m_contact_cols::m_p1_id) =
        std::get<m_contact_cols::m_p1_id>(elem);
    m_c_data(counter, m_contact_cols::m_p2_id) =
        std::get<m_contact_cols::m_p2_id>(elem);
    m_c_data(counter, m_contact_cols::m_contact_overlap) =
        std::get<m_contact_cols::m_contact_overlap>(elem);
    m_c_data(counter, m_contact_cols::m_cn_force_x) =
        std::get<m_contact_cols::m_cn_force_x>(elem);
    m_c_data(counter, m_contact_cols::m_cn_force_y) =
        std::get<m_contact_cols::m_cn_force_y>(elem);
    m_c_data(counter, m_contact_cols::m_cn_force_z) =
        std::get<m_contact_cols::m_cn_force_z>(elem);
    m_c_data(counter, m_contact_cols::m_ct_force_x) =
        std::get<m_contact_cols::m_ct_force_x>(elem);
    m_c_data(counter, m_contact_cols::m_ct_force_y) =
        std::get<m_contact_cols::m_ct_force_y>(elem);
    m_c_data(counter, m_contact_cols::m_ct_force_z) =
        std::get<m_contact_cols::m_ct_force_z>(elem);
    m_c_data(counter, m_contact_cols::m_c_force_x) =
        std::get<m_contact_cols::m_c_force_x>(elem);
    m_c_data(counter, m_contact_cols::m_c_force_y) =
        std::get<m_contact_cols::m_c_force_y>(elem);
    m_c_data(counter, m_contact_cols::m_c_force_z) =
        std::get<m_contact_cols::m_c_force_z>(elem);
    m_c_data(counter, m_contact_cols::m_p1_r) =
        radius.at(std::get<m_contact_cols::m_p1_id>(elem));
    m_c_data(counter, m_contact_cols::m_p2_r) =
        radius.at(std::get<m_contact_cols::m_p2_id>(elem));
    m_c_data(counter, m_contact_cols::m_c_sliding_contact) =
        std::get<m_contact_cols::m_c_sliding_contact>(elem);
    calc_pot_energy(counter);
    counter++;
  }
  counter = 0;
  m_p_data = d_mat(vec_p_size, m_particle_cols_items);
  for (const auto &elem : particle_data) {
    p_cellstridx.insert(std::make_pair(std::get<13>(elem), counter));
    m_p_data(counter, m_part_cols::m_p_id) =
        std::get<m_part_cols::m_p_id>(elem);
    m_p_data(counter, m_part_cols::m_p_vx) =
        std::get<m_part_cols::m_p_vx>(elem);
    m_p_data(counter, m_part_cols::m_p_vy) =
        std::get<m_part_cols::m_p_vy>(elem);
    m_p_data(counter, m_part_cols::m_p_vz) =
        std::get<m_part_cols::m_p_vz>(elem);
    m_p_data(counter, m_part_cols::m_p_coord) =
        std::get<m_part_cols::m_p_coord>(elem);
    m_p_data(counter, m_part_cols::m_p_disp_x) =
        std::get<m_part_cols::m_p_disp_x>(elem);
    m_p_data(counter, m_part_cols::m_p_disp_y) =
        std::get<m_part_cols::m_p_disp_y>(elem);
    m_p_data(counter, m_part_cols::m_p_disp_z) =
        std::get<m_part_cols::m_p_disp_z>(elem);
    m_p_data(counter, m_part_cols::m_p_ke_rot) =
        std::get<m_part_cols::m_p_ke_rot>(elem);
    m_p_data(counter, m_part_cols::m_p_ke_tra) =
        std::get<m_part_cols::m_p_ke_tra>(elem);
    m_p_data(counter, m_part_cols::m_p_omegax) =
        std::get<m_part_cols::m_p_omegax>(elem);
    m_p_data(counter, m_part_cols::m_p_omegay) =
        std::get<m_part_cols::m_p_omegay>(elem);
    m_p_data(counter, m_part_cols::m_p_omegaz) =
        std::get<m_part_cols::m_p_omegaz>(elem);
    counter++;
  }
  BOOST_LOG_TRIVIAL(info) << "Initalized Matrix ts: " << ts;
}

inline void PropertyCalculator::calc_pot_energy(int &counter)
{
  m_c_data(counter, m_contact_cols::m_radstar) =
      m_c_data(counter, m_contact_cols::m_p1_r) *
      m_c_data(counter, m_contact_cols::m_p2_r) /
      (m_c_data(counter, m_contact_cols::m_p1_r) +
       m_c_data(counter, m_contact_cols::m_p2_r));
  if (m_c_data(counter, m_contact_cols::m_p1_r) <= 0)
    m_c_data(counter, m_contact_cols::m_radstar) =
        m_c_data(counter, m_contact_cols::m_p2_r);
  if (m_c_data(counter, m_contact_cols::m_p2_r) <= 0)
    m_c_data(counter, m_contact_cols::m_radstar) =
        m_c_data(counter, m_contact_cols::m_p1_r);
  m_c_data(counter, m_contact_cols::m_kn) =
      4 / 3.0 * Ystar *
      std::sqrt(m_c_data(counter, m_contact_cols::m_radstar) *
                m_c_data(counter, m_contact_cols::m_contact_overlap));
  m_c_data(counter, m_contact_cols::m_kt) =
      8.0 * Gstar *
      std::sqrt(m_c_data(counter, m_contact_cols::m_radstar) *
                m_c_data(counter, m_contact_cols::m_contact_overlap));
  m_c_data(counter, m_contact_cols::m_fnor) =
      std::sqrt(std::pow(m_c_data(counter, m_contact_cols::m_cn_force_x), 2) +
                std::pow(m_c_data(counter, m_contact_cols::m_cn_force_y), 2) +
                std::pow(m_c_data(counter, m_contact_cols::m_cn_force_z), 2));
  m_c_data(counter, m_contact_cols::m_ftan) =
      std::sqrt(std::pow(m_c_data(counter, m_contact_cols::m_ct_force_x), 2) +
                std::pow(m_c_data(counter, m_contact_cols::m_ct_force_y), 2) +
                std::pow(m_c_data(counter, m_contact_cols::m_ct_force_z), 2));
  m_c_data(counter, m_contact_cols::m_penor) =
      0.5 * (std::pow(m_c_data(counter, m_contact_cols::m_fnor), 2) /
             m_c_data(counter, m_contact_cols::m_kn));
  m_c_data(counter, m_contact_cols::m_petan) =
      0.5 * (std::pow(m_c_data(counter, m_contact_cols::m_ftan), 2) /
             m_c_data(counter, m_contact_cols::m_kt));
  if (m_c_data(counter, m_contact_cols::m_contact_overlap) >=
      0.1 * m_c_data(counter, m_contact_cols::m_radstar)) {
    m_c_data(counter, m_contact_cols::m_fnor) = 0.0;
    m_c_data(counter, m_contact_cols::m_ftan) = 0.0;
    m_c_data(counter, m_contact_cols::m_penor) = 0.0;
    m_c_data(counter, m_contact_cols::m_petan) = 0.0;
  }
}

void PropertyCalculator::calc_global_state()
{
  auto c_rows = m_c_data.rows();
  global_map["penor"] = m_c_data.col(m_contact_cols::m_penor).sum();
  global_map["petan"] = m_c_data.col(m_contact_cols::m_petan).sum();
  global_map["ftan"] = m_c_data.col(m_contact_cols::m_ftan).sum();
  global_map["fnor"] = m_c_data.col(m_contact_cols::m_fnor).sum();
  global_map["slipping_ration"] =
      m_c_data.col(m_contact_cols::m_c_sliding_contact).sum() / c_rows;
  global_map["cforce_x"] =
      m_c_data.col(m_contact_cols::m_c_force_x).sum() / c_rows;
  global_map["cforce_y"] =
      m_c_data.col(m_contact_cols::m_c_force_y).sum() / c_rows;
  global_map["cforce_z"] =
      m_c_data.col(m_contact_cols::m_c_force_z).sum() / c_rows;

  global_map["cforce_mag"] =
      Eigen::sqrt(
          Eigen::square(m_c_data.col(m_contact_cols::m_c_force_x).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_c_force_y).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_c_force_z).array()))
          .sum() /
      c_rows;
  m_c_data.col(m_contact_cols::m_c_force_mag).sum();
  global_map["ct_force_x"] =
      m_c_data.col(m_contact_cols::m_ct_force_x).sum() / c_rows;
  global_map["ct_force_y"] =
      m_c_data.col(m_contact_cols::m_ct_force_y).sum() / c_rows;
  global_map["ct_force_z"] =
      m_c_data.col(m_contact_cols::m_ct_force_z).sum() / c_rows;
  global_map["ct_force_mag"] =
      Eigen::sqrt(
          Eigen::square(m_c_data.col(m_contact_cols::m_ct_force_x).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_ct_force_y).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_ct_force_z).array()))
          .sum() /
      c_rows;
  global_map["cn_force_x"] =
      m_c_data.col(m_contact_cols::m_cn_force_x).sum() / c_rows;
  global_map["cn_force_y"] =
      m_c_data.col(m_contact_cols::m_cn_force_y).sum() / c_rows;
  global_map["cn_force_z"] =
      m_c_data.col(m_contact_cols::m_cn_force_z).sum() / c_rows;
  global_map["cn_force_mag"] =
      Eigen::sqrt(
          Eigen::square(m_c_data.col(m_contact_cols::m_cn_force_x).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_cn_force_y).array()) +
          Eigen::square(m_c_data.col(m_contact_cols::m_cn_force_z).array()))
          .sum() /
      c_rows;
  auto p_rows = m_p_data.rows();
  global_map["vel_x"] = m_p_data.col(m_part_cols::m_p_vx).sum() / p_rows;
  global_map["vel_y"] = m_p_data.col(m_part_cols::m_p_vy).sum() / p_rows;
  global_map["vel_z"] = m_p_data.col(m_part_cols::m_p_vz).sum() / p_rows;
  global_map["vel_mag"] =
      Eigen::sqrt(Eigen::square(m_p_data.col(m_part_cols::m_p_vx).array()) +
                  Eigen::square(m_p_data.col(m_part_cols::m_p_vy).array()) +
                  Eigen::square(m_p_data.col(m_part_cols::m_p_vz).array()))
          .sum() /
      p_rows;
  global_map["coord"] = m_p_data.col(m_part_cols::m_p_coord).sum() / p_rows;
  global_map["ke_rot"] = m_p_data.col(m_part_cols::m_p_ke_rot).sum();
  global_map["ke_tra"] = m_p_data.col(m_part_cols::m_p_ke_tra).sum();
  global_map["ke_tot"] = global_map["ke_rot"] + global_map["ke_tra"];
  global_map["disp_x"] = m_p_data.col(m_part_cols::m_p_disp_x).sum() / p_rows;
  global_map["disp_y"] = m_p_data.col(m_part_cols::m_p_disp_y).sum() / p_rows;
  global_map["disp_z"] = m_p_data.col(m_part_cols::m_p_disp_z).sum() / p_rows;
  global_map["disp_mag"] =
      Eigen::sqrt(Eigen::square(m_p_data.col(m_part_cols::m_p_disp_x).array()) +
                  Eigen::square(m_p_data.col(m_part_cols::m_p_disp_y).array()) +
                  Eigen::square(m_p_data.col(m_part_cols::m_p_disp_z).array()))
          .sum() /
      p_rows;
  for (auto &elem : global_map) {
    BOOST_LOG_TRIVIAL(info) << elem.first << ": " << elem.second;
  }
}

void PropertyCalculator::aggregate_per_cell()
{
  for (auto &elem : decomp_str) {
    auto c_range = c_cellstridx.equal_range(elem.cellstr);
    Accumulator<double> acc_penor, acc_petan, acc_ftan, acc_fnor, acc_slipping,
        acc_cforce_x, acc_cforce_y, acc_cforce_z, acc_ct_force_x,
        acc_ct_force_y, acc_ct_force_z, acc_cn_force_x, acc_cn_force_y,
        acc_cn_force_z, acc_c_force_mag, acc_ct_force_mag, acc_cn_force_mag;
    if (c_range.first != c_range.second) {
      for (auto it = c_range.first; it != c_range.second; ++it) {
        acc_penor(m_c_data((*it).second, m_contact_cols::m_penor));
        acc_petan(m_c_data((*it).second, m_contact_cols::m_petan));
        acc_ftan(m_c_data((*it).second, m_contact_cols::m_ftan));
        acc_fnor(m_c_data((*it).second, m_contact_cols::m_fnor));
        acc_slipping(
            m_c_data((*it).second, m_contact_cols::m_c_sliding_contact));
        acc_cforce_x(m_c_data((*it).second, m_contact_cols::m_c_force_x));
        acc_cforce_y(m_c_data((*it).second, m_contact_cols::m_c_force_y));
        acc_cforce_z(m_c_data((*it).second, m_contact_cols::m_c_force_z));
        acc_ct_force_x(m_c_data((*it).second, m_contact_cols::m_ct_force_x));
        acc_ct_force_y(m_c_data((*it).second, m_contact_cols::m_ct_force_y));
        acc_ct_force_z(m_c_data((*it).second, m_contact_cols::m_ct_force_z));
        acc_cn_force_x(m_c_data((*it).second, m_contact_cols::m_cn_force_x));
        acc_cn_force_y(m_c_data((*it).second, m_contact_cols::m_cn_force_y));
        acc_cn_force_z(m_c_data((*it).second, m_contact_cols::m_cn_force_z));
        acc_c_force_mag(std::sqrt(
            std::pow(m_c_data((*it).second, m_contact_cols::m_c_force_x), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_c_force_y), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_c_force_z), 2)));
        acc_ct_force_mag(std::sqrt(
            std::pow(m_c_data((*it).second, m_contact_cols::m_ct_force_x), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_ct_force_y), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_ct_force_z), 2)));
        acc_cn_force_mag(std::sqrt(
            std::pow(m_c_data((*it).second, m_contact_cols::m_cn_force_x), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_cn_force_y), 2) +
            std::pow(m_c_data((*it).second, m_contact_cols::m_cn_force_z), 2)));
      }
    }
    auto p_range = p_cellstridx.equal_range(elem.cellstr);
    Accumulator<double> acc_vel_mag, acc_vel_x, acc_vel_y, acc_vel_z, acc_coord,
        acc_ke_rot, acc_ke_tra, acc_ke_tot, acc_disp_x, acc_disp_y, acc_disp_z,
        acc_disp_mag;
    if (p_range.first != p_range.second) {
      for (auto it = p_range.first; it != p_range.second; ++it) {
        acc_vel_x(m_p_data((*it).second, m_part_cols::m_p_vx));
        acc_vel_y(m_p_data((*it).second, m_part_cols::m_p_vy));
        acc_vel_z(m_p_data((*it).second, m_part_cols::m_p_vz));
        acc_vel_mag(std::sqrt(
            std::pow(m_p_data((*it).second, m_part_cols::m_p_vx), 2) +
            std::pow(m_p_data((*it).second, m_part_cols::m_p_vy), 2) +
            std::pow(m_p_data((*it).second, m_part_cols::m_p_vz), 2)));
        acc_coord(m_p_data((*it).second, m_part_cols::m_p_coord));
        acc_ke_rot(m_p_data((*it).second, m_part_cols::m_p_ke_rot));
        acc_ke_tra(m_p_data((*it).second, m_part_cols::m_p_ke_tra));
        acc_ke_tot(m_p_data((*it).second, m_part_cols::m_p_ke_rot) +
                   m_p_data((*it).second, m_part_cols::m_p_ke_tra));
        acc_disp_x(m_p_data((*it).second, m_part_cols::m_p_disp_x));
        acc_disp_y(m_p_data((*it).second, m_part_cols::m_p_disp_y));
        acc_disp_z(m_p_data((*it).second, m_part_cols::m_p_disp_z));
        acc_disp_mag(std::sqrt(
            std::pow(m_p_data((*it).second, m_part_cols::m_p_disp_x), 2) +
            std::pow(m_p_data((*it).second, m_part_cols::m_p_disp_y), 2) +
            std::pow(m_p_data((*it).second, m_part_cols::m_p_disp_z), 2)));
      }
    }

    std::map<std::string, double> agg_elem{};
    agg_elem["penor"] = acc_penor.getAccVal();
    agg_elem["petan"] = acc_petan.getAccVal();
    agg_elem["ftan"] = acc_ftan.getAccVal();
    agg_elem["fnor"] = acc_fnor.getAccVal();
    agg_elem["slipping_ratio"] = acc_slipping.get_mean();
    agg_elem["cforce_x"] = acc_cforce_x.get_mean();
    agg_elem["cforce_y"] = acc_cforce_y.get_mean();
    agg_elem["cforce_z"] = acc_cforce_z.get_mean();
    agg_elem["ct_force_x"] = acc_ct_force_x.get_mean();
    agg_elem["ct_force_y"] = acc_ct_force_y.get_mean();
    agg_elem["ct_force_z"] = acc_ct_force_z.get_mean();
    agg_elem["cn_force_x"] = acc_cn_force_x.get_mean();
    agg_elem["cn_force_y"] = acc_cn_force_y.get_mean();
    agg_elem["cn_force_z"] = acc_cn_force_z.get_mean();
    agg_elem["c_force_mag"] = acc_c_force_mag.get_mean();
    agg_elem["ct_force_mag"] = acc_ct_force_mag.get_mean();
    agg_elem["cn_force_mag"] = acc_cn_force_mag.get_mean();
    agg_elem["vel_x"] = acc_vel_x.get_mean();
    agg_elem["vel_y"] = acc_vel_y.get_mean();
    agg_elem["vel_z"] = acc_vel_z.get_mean();
    agg_elem["vel_mag"] = acc_vel_mag.get_mean();
    agg_elem["coord"] = acc_coord.get_mean();
    agg_elem["ke_rot"] = acc_ke_rot.getAccVal();
    agg_elem["ke_tra"] = acc_ke_tra.getAccVal();
    agg_elem["ke_tot"] = acc_ke_tot.getAccVal();
    agg_elem["disp_x"] = acc_disp_x.get_mean();
    agg_elem["disp_y"] = acc_disp_y.get_mean();
    agg_elem["disp_z"] = acc_disp_z.get_mean();
    agg_elem["disp_mag"] = acc_disp_mag.get_mean();

    aggregate_map[elem.cellstr] = agg_elem;
  }

  calc_global_state();
}
const aggregate_map_t &PropertyCalculator::getAggregateMap() const
{
  return aggregate_map;
}

PropertyCalculator::~PropertyCalculator() {}

const std::map<std::string, double> &PropertyCalculator::getGlobalMap() const
{
  return global_map;
}
