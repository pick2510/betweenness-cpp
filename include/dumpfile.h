#ifndef DUMPFILE_H
#define DUMPFILE_H

#include <fstream>
#include <string>
#include <vector>
#include "data.h"

class dumpfile {
  static constexpr int begin_line = 9;
  static constexpr int ts_line = 1;

private:
  std::ifstream file;
  std::vector<ContactColumns> ts_file_column{};
  int timestep;
  void set_fpointer(int n);
public:
  dumpfile(const std::string &Path);
  void parse_file();
  std::vector<ContactColumns> getData(){return ts_file_column;}
  ~dumpfile() { file.close(); };
};

#endif