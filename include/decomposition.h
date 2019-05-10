#ifndef DECOMP_H
#define DECOMP_H

#include "data.h"
#include <set>
#include <string>
#include <vector>

class Decomposition {
private:
  const Config runningConf;
  double dx, dy, dz;
  std::vector<std::string> cells;
  std::vector<int> x_range, y_range, z_range;

public:
  Decomposition(const Config &runningConf);
  std::string calc_cell(coordinate &coord) const;
  cell calc_cell_numeric(coordinate &coord) const;
  ~Decomposition();
  const std::vector<std::string> &getCells() const;
};

#endif