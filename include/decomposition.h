#ifndef DECOMP_H
#define DECOMP_H

#include "data.h"
#include <string>
#include <vector>
#include <set>

class Decomposition {
private:
  const Config runningConf;
  double dx, dy, dz;
  std::vector<cell> cells;
  std::vector<int> x_range,y_range, z_range;

public:
  Decomposition(const Config &runningConf);
  std::string calc_cell(double x, double y, double z);
  cell calc_cell_numeric(double x, double y, double z);
  ~Decomposition();
};

#endif