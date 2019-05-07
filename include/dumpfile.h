#ifndef DUMPFILE_H
#define DUMPFILE_H

#include "data.h"
#include "decomposition.h"
#include "utils.h"
#include <fstream>
#include <string>
#include <vector>

class dumpfile {
  static constexpr int begin_line = 9;
  static constexpr int ts_line = 1;

private:
  std::ifstream file;
  dumpfile_t type;
  const std::map<int, double> radius;
  const Decomposition decomp;
  std::vector<ContactColumns> ts_file_column{};
  std::vector<ParticleColumns> part_file_column{};
  long timestep;
  inline void set_fpointer(int n) { goto_line(file, n); }

public:
  dumpfile(const std::string &Path, const std::map<int, double> &radius,
           const Decomposition &decomp, dumpfile_t type);
  void parse_file();
  void parse_contacts();
  void parse_particles();
  long gettimestep() const { return timestep; }
  inline std::vector<ContactColumns> getContactData() { return ts_file_column; }
  inline std::vector<ParticleColumns> getParticleData()
  {
    return part_file_column;
  }
  static coordinate calc_contactpoint(ContactColumns &contact,
                                      const std::map<int, double> &radius);

  ~dumpfile() { file.close(); };
};

#endif