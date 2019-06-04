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
    : radius(radius), contact_data(contact_data), vec_size(contact_data.size()),
      decomp_str(decomp_str)
{
  tuple_size = std::tuple_size<std::decay_t<decltype(contact_data[0])>>::value;
  m_data = d_mat(vec_size, m_contact_cols_items);
  ts = std::get<10>(contact_data[0]);
  int counter = 0;
  for (const auto &elem : contact_data) {
    cellstridx.insert(std::pair<std::string, int>(std::get<9>(elem), counter));
    m_data(counter, m_contact_cols::m_p1_id) =
        std::get<m_contact_cols::m_p1_id>(elem);
    m_data(counter, m_contact_cols::m_p2_id) =
        std::get<m_contact_cols::m_p2_id>(elem);
    m_data(counter, m_contact_cols::m_contact_overlap) =
        std::get<m_contact_cols::m_contact_overlap>(elem);
    m_data(counter, m_contact_cols::m_cn_force_x) =
        std::get<m_contact_cols::m_cn_force_x>(elem);
    m_data(counter, m_contact_cols::m_cn_force_y) =
        std::get<m_contact_cols::m_cn_force_y>(elem);
    m_data(counter, m_contact_cols::m_cn_force_z) =
        std::get<m_contact_cols::m_cn_force_z>(elem);
    m_data(counter, m_contact_cols::m_ct_force_x) =
        std::get<m_contact_cols::m_ct_force_x>(elem);
    m_data(counter, m_contact_cols::m_ct_force_y) =
        std::get<m_contact_cols::m_ct_force_y>(elem);
    m_data(counter, m_contact_cols::m_ct_force_z) =
        std::get<m_contact_cols::m_ct_force_z>(elem);
    m_data(counter, m_contact_cols::m_p1_r) =
        radius.at(std::get<m_contact_cols::m_p1_id>(elem));
    m_data(counter, m_contact_cols::m_p2_r) =
        radius.at(std::get<m_contact_cols::m_p2_id>(elem));
    calc_pot_energy(counter);
    counter++;
  }
  BOOST_LOG_TRIVIAL(info) << "Initalized Matrix ts: " << ts;
}

inline void PropertyCalculator::calc_pot_energy(int &counter)
{
  m_data(counter, m_contact_cols::m_radstar) =
      m_data(counter, m_contact_cols::m_p1_r) *
      m_data(counter, m_contact_cols::m_p2_r) /
      (m_data(counter, m_contact_cols::m_p1_r) +
       m_data(counter, m_contact_cols::m_p2_r));
  if (m_data(counter, m_contact_cols::m_p1_r) <= 0)
    m_data(counter, m_contact_cols::m_radstar) =
        m_data(counter, m_contact_cols::m_p2_r);
  if (m_data(counter, m_contact_cols::m_p2_r) <= 0)
    m_data(counter, m_contact_cols::m_radstar) =
        m_data(counter, m_contact_cols::m_p1_r);
  m_data(counter, m_contact_cols::m_kn) =
      4 / 3.0 * Ystar *
      std::sqrt(m_data(counter, m_contact_cols::m_radstar) *
                m_data(counter, m_contact_cols::m_contact_overlap));
  m_data(counter, m_contact_cols::m_kt) =
      8.0 * Gstar *
      std::sqrt(m_data(counter, m_contact_cols::m_radstar) *
                m_data(counter, m_contact_cols::m_contact_overlap));
  m_data(counter, m_contact_cols::m_fnor) =
      std::sqrt(std::pow(m_data(counter, m_contact_cols::m_cn_force_x), 2) +
                std::pow(m_data(counter, m_contact_cols::m_cn_force_y), 2) +
                std::pow(m_data(counter, m_contact_cols::m_cn_force_z), 2));
  m_data(counter, m_contact_cols::m_ftan) =
      std::sqrt(std::pow(m_data(counter, m_contact_cols::m_ct_force_x), 2) +
                std::pow(m_data(counter, m_contact_cols::m_ct_force_y), 2) +
                std::pow(m_data(counter, m_contact_cols::m_ct_force_z), 2));
  m_data(counter, m_contact_cols::m_penor) =
      0.5 * (std::pow(m_data(counter, m_contact_cols::m_fnor), 2) /
             m_data(counter, m_contact_cols::m_kn));
  m_data(counter, m_contact_cols::m_petan) =
      0.5 * (std::pow(m_data(counter, m_contact_cols::m_ftan), 2) /
             m_data(counter, m_contact_cols::m_kt));
  if (m_data(counter, m_contact_cols::m_contact_overlap) >=
      0.1 * m_data(counter, m_contact_cols::m_radstar)) {
    m_data(counter, m_contact_cols::m_fnor) = 0.0;
    m_data(counter, m_contact_cols::m_ftan) = 0.0;
    m_data(counter, m_contact_cols::m_penor) = 0.0;
    m_data(counter, m_contact_cols::m_petan) = 0.0;
  }
}

void PropertyCalculator::aggregate_per_cell()
{
  for (auto &elem : decomp_str) {
    auto range = cellstridx.equal_range(elem.cellstr);
    if (range.first == range.second)
      continue;
    Accumulator<double> acc_penor, acc_petan, acc_ftan, acc_fnor;
    for (auto it = range.first; it != range.second; ++it) {
      acc_penor(m_data((*it).second, m_contact_cols::m_penor));
      acc_petan(m_data((*it).second, m_contact_cols::m_petan));
      acc_ftan(m_data((*it).second, m_contact_cols::m_ftan));
      acc_fnor(m_data((*it).second, m_contact_cols::m_fnor));
    }
    std::map<std::string, double> agg_elem{};
    agg_elem["penor"] = acc_penor.getAccVal();
    agg_elem["petan"] = acc_petan.getAccVal();
    agg_elem["ftan"] = acc_ftan.getAccVal();
    agg_elem["fnor"] = acc_fnor.getAccVal();
    aggregate_map[elem.cellstr] = agg_elem;
  }
}
const aggregate_map_t &PropertyCalculator::getAggregateMap() const
{
  return aggregate_map;
}

PropertyCalculator::~PropertyCalculator() {}
