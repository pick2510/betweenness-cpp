#ifndef DUMPFILE_H
#define DUMPFILE_H

#include <fstream>
#include <string>
#include <vector>
#include "data.h"
#include "utils.h"
#include "decomposition.h"

class dumpfile {
  static constexpr int begin_line = 9;
  static constexpr int ts_line = 1;

private:
  std::ifstream file;
  const std::vector<double> radius;
  const Decomposition decomp;
  std::vector<ContactColumns> ts_file_column{};
  long timestep;
  inline void set_fpointer(int n) { goto_line(file, n); }
  coordinate calc_contactpoint(ContactColumns &contact);

public:
  dumpfile(const std::string &Path, const std::vector<double> &radius, const Decomposition &decomp);
  void parse_file();
  long gettimestep() const {return timestep;}
  inline std::vector<ContactColumns> getData(){return ts_file_column;}
  ~dumpfile() { file.close(); };
};

#endif