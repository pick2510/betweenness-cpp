#include "dumpfile.h"
#include "data.h"
#include "decomposition.h"
#include "sqlite_orm.h"
#include "utils.h"
#include <boost/log/trivial.hpp>

dumpfile::dumpfile(const std::string &Path, const std::vector<double> &radius, const Decomposition &decomp)
    : file(Path, std::ios::in), radius(radius), decomp(decomp)
{
  set_fpointer(dumpfile::ts_line);
  file >> timestep;
  BOOST_LOG_TRIVIAL(info) << timestep;
  set_fpointer(dumpfile::begin_line);
  BOOST_LOG_TRIVIAL(info) << "Initialized";
}



void dumpfile::parse_file()
{
  std::string line;
  while (std::getline(file, line)) {
    ContactColumns columns{};
    std::vector<std::string> splitted_line;
    split_string(line, splitted_line);
    columns.p1_x = std::stod(splitted_line[ContactTXTColumns::p1_x]);
    columns.p1_y = std::stod(splitted_line[ContactTXTColumns::p1_y]);
    columns.p1_z = std::stod(splitted_line[ContactTXTColumns::p1_z]);
    columns.p2_x = std::stod(splitted_line[ContactTXTColumns::p2_x]);
    columns.p2_y = std::stod(splitted_line[ContactTXTColumns::p2_y]);
    columns.p2_z = std::stod(splitted_line[ContactTXTColumns::p2_z]);
    columns.p1_vx = std::stod(splitted_line[ContactTXTColumns::p1_vx]);
    columns.p1_vy = std::stod(splitted_line[ContactTXTColumns::p1_vy]);
    columns.p1_vz = std::stod(splitted_line[ContactTXTColumns::p1_vz]);
    columns.p2_vx = std::stod(splitted_line[ContactTXTColumns::p2_vx]);
    columns.p2_vy = std::stod(splitted_line[ContactTXTColumns::p2_vy]);
    columns.p2_vz = std::stod(splitted_line[ContactTXTColumns::p2_vz]);
    columns.p1_id = std::stoi(splitted_line[ContactTXTColumns::p1_id]);
    columns.p2_id = std::stoi(splitted_line[ContactTXTColumns::p2_id]);
    columns.is_periodic =
        std::stoi(splitted_line[ContactTXTColumns::is_periodic]);
    columns.c_force_x = std::stod(splitted_line[ContactTXTColumns::c_force_x]);
    columns.c_force_y = std::stod(splitted_line[ContactTXTColumns::c_force_y]);
    columns.c_force_z = std::stod(splitted_line[ContactTXTColumns::c_force_z]);
    columns.cn_force_x =
        std::stod(splitted_line[ContactTXTColumns::cn_force_x]);
    columns.cn_force_y =
        std::stod(splitted_line[ContactTXTColumns::cn_force_y]);
    columns.cn_force_z =
        std::stod(splitted_line[ContactTXTColumns::cn_force_z]);
    columns.ct_force_x =
        std::stod(splitted_line[ContactTXTColumns::ct_force_x]);
    columns.ct_force_y =
        std::stod(splitted_line[ContactTXTColumns::ct_force_y]);
    columns.ct_force_z =
        std::stod(splitted_line[ContactTXTColumns::ct_force_z]);
    columns.c_torque_x =
        std::stod(splitted_line[ContactTXTColumns::c_torque_x]);
    columns.c_torque_y =
        std::stod(splitted_line[ContactTXTColumns::c_torque_y]);
    columns.c_torque_z =
        std::stod(splitted_line[ContactTXTColumns::c_torque_z]);
    columns.disp_x = std::stod(splitted_line[ContactTXTColumns::disp_x]);
    columns.disp_y = std::stod(splitted_line[ContactTXTColumns::disp_y]);
    columns.disp_z = std::stod(splitted_line[ContactTXTColumns::disp_z]);
    columns.contact_area =
        std::stod(splitted_line[ContactTXTColumns::contact_area]);
    columns.contact_overlap =
        std::stod(splitted_line[ContactTXTColumns::contact_overlap]);
    columns.sliding_contact =
        std::stod(splitted_line[ContactTXTColumns::sliding_contact]);
    columns.ts = timestep;
    auto coord = calc_contactpoint(columns);
    auto cell = decomp.calc_cell_numeric(coord);
    columns.cellstr = decomp.calc_cell(coord);
    columns.cell_x = cell.x;
    columns.cell_y = cell.y;
    columns.cell_z = cell.z;
    ts_file_column.push_back(columns);
  }
}

coordinate dumpfile::calc_contactpoint(ContactColumns &contact){
    double Ab_x, Ab_y, Ab_z, cm_x, cm_y, cm_z;
    double p1_r {radius[contact.p1_id]};
    double p2_r {radius[contact.p2_id]};
    if (p1_r < 0) p1_r = p2_r;
    if (p2_r < 0) p2_r = p1_r;
    double da = p1_r / (p1_r + p2_r);
    Ab_x = contact.p2_x - contact.p1_x;
    Ab_y = contact.p2_y - contact.p1_y;
    Ab_z = contact.p2_z - contact.p1_z;
    cm_x = p1_x + Ab_x * da;
    cm_y = p1_y + Ab_y * da;
    cm_z = p1_z + Ab_z * da;
    return coordinate{
        .x = cm_x,
        .y = cm_y,
        .z = cm_z
    };
}