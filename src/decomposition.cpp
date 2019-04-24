#include "decomposition.h"
#include "data.h"
#include "utils.h"
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>

Decomposition::Decomposition(const Config &runningConf)
    : runningConf(runningConf), x_range(runningConf.x_cells, 0),
      y_range(runningConf.y_cells, 0), z_range(runningConf.z_cells, 0)
{

  dx = runningConf.domainsize_x / runningConf.x_cells;
  dy = runningConf.domainsize_y / runningConf.y_cells;
  dz = runningConf.domainsize_z / runningConf.z_cells;
  std::iota(x_range.begin(), x_range.end(), 0);
  std::iota(y_range.begin(), y_range.end(), 0);
  std::iota(z_range.begin(), z_range.end(), 0);
  cells = getCartesianProduct(x_range, y_range, z_range);
}

std::string Decomposition::calc_cell(double x, double y, double z)
{
  std::stringstream ss;
  ss << "cell" << std::setfill('0') << std::setw(2) << std::floor(x / dx) << "_"
     << std::setfill('0') << std::setw(2) << std::floor(y / dy) << "_"
     << std::setfill('0') << std::setw(2) << std::floor(z / dz);
  return ss.str();
}

cell Decomposition::calc_cell_numeric(double x, double y, double z){
    return cell{
        .x = static_cast<int>(std::floor(x / dx)),
        .y = static_cast<int>(std::floor(y / dy)),
        .z = static_cast<int>(std::floor(z / dz))
    };
}


Decomposition::~Decomposition() {}
