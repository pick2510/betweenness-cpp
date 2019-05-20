#ifndef DATA_H
#define DATA_H
#include <algorithm>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <cmath>
#include <map>
#include <string>
#include <vector>

constexpr int MASTER = 0;
constexpr int TAG_RESULT = 1;
constexpr int TAG_BREAK = 2;
constexpr int TAG_TS = 5;
constexpr int TAG_SIZE = 3;
constexpr int TAG_FILE = 10;
constexpr int HOSTNAME_LEN = 255;
constexpr int TAG_PART_TS = 20;

struct cell {
  int x, y, z;
};

struct coordinate {
  double x, y, z;
};

struct Config {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &InputPath;
    ar &OutputPath;
    ar &sep;
    ar &x_cells;
    ar &y_cells;
    ar &z_cells;
    ar &domainsize_x;
    ar &domainsize_y;
    ar &domainsize_z;
    ar &output_percentile;
    ar &randomly_selected;
    ar &contact_filename;
    ar &particle_filename;
    ar &v_particles;
    ar &spinup_time;
    ar &chunk_len;
  }

public:
  std::string InputPath;
  std::string OutputPath;
  char sep[2] = " ";
  int x_cells = 0;
  int y_cells = 0;
  int z_cells = 0;
  double domainsize_x = 0.0;
  double domainsize_y = 0.0;
  double domainsize_z = 0.0;
  double output_percentile = 0.9;
  int randomly_selected = 20;
  long spinup_time;
  int chunk_len;
  std::string contact_filename;
  std::string particle_filename;
  std::string v_particles;
};

struct Result {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &ts;
    ar &keys;
    ar &vals;
    ar &mean;
    ar &var;
    ar &std;
    ar &skew;
    ar &kur;
    ar &q_090;
    ar &q_099;
  }

public:
  long ts;
  std::vector<int> keys;
  std::vector<double> vals;
  double mean;
  double var;
  double std;
  double skew;
  double kur;
  double q_090;
  double q_099;
};

constexpr char ts_particle_path[] = "particles";
constexpr char ts_centrality_path[] = "centrality";

struct ContactColumns {
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &id;
    ar &p1_x;
    ar &p1_y;
    ar &p1_z;
    ar &p2_x;
    ar &p2_y;
    ar &p2_z;
    ar &p1_vx;
    ar &p1_vy;
    ar &p1_vz;
    ar &p2_vx;
    ar &p2_vy;
    ar &p2_vz;
    ar &p1_id;
    ar &p2_id;
    ar &is_periodic;
    ar &c_force_x;
    ar &c_force_y;
    ar &c_force_z;
    ar &cn_force_x;
    ar &cn_force_y;
    ar &cn_force_z;
    ar &ct_force_x;
    ar &ct_force_y;
    ar &ct_force_z;
    ar &c_torque_x;
    ar &c_torque_y;
    ar &c_torque_z;
    ar &disp_x;
    ar &disp_y;
    ar &disp_z;
    ar &contact_area;
    ar &contact_overlap;
    ar &sliding_contact;
    ar &cell_x;
    ar &cell_y;
    ar &cell_z;
    ar &ts;
    ar &cellstr;
  }

  long id;
  double p1_x;
  double p1_y;
  double p1_z;
  double p2_x;
  double p2_y;
  double p2_z;
  double p1_vx;
  double p1_vy;
  double p1_vz;
  double p2_vx;
  double p2_vy;
  double p2_vz;
  int p1_id;
  int p2_id;
  int is_periodic;
  double c_force_x;
  double c_force_y;
  double c_force_z;
  double cn_force_x;
  double cn_force_y;
  double cn_force_z;
  double ct_force_x;
  double ct_force_y;
  double ct_force_z;
  double c_torque_x;
  double c_torque_y;
  double c_torque_z;
  double disp_x;
  double disp_y;
  double disp_z;
  double contact_area;
  double contact_overlap;
  int sliding_contact;
  int cell_x;
  int cell_y;
  int cell_z;
  long ts;
  std::string cellstr;
};

struct ParticleColumns {
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &id;
    ar &p_id;
    ar &p_type;
    ar &p_x;
    ar &p_y;
    ar &p_z;
    ar &p_rad;
    ar &p_vx;
    ar &p_vy;
    ar &p_vz;
    ar &p_fx;
    ar &p_fy;
    ar &p_fz;
    ar &p_omegax;
    ar &p_omegay;
    ar &p_omegaz;
    ar &p_coord;
    ar &p_disp_x;
    ar &p_disp_y;
    ar &p_disp_z;
    ar &p_disp_mag;
    ar &p_ke_rot;
    ar &p_ke_tra;
    ar &cell_x;
    ar &cell_y;
    ar &cell_z;
    ar &ts;
    ar &cellstr;
  }

  long id;
  int p_id;
  int p_type;
  double p_x;
  double p_y;
  double p_z;
  double p_rad;
  double p_vx;
  double p_vy;
  double p_vz;
  double p_fx;
  double p_fy;
  double p_fz;
  double p_omegax;
  double p_omegay;
  double p_omegaz;
  int p_coord;
  double p_disp_x;
  double p_disp_y;
  double p_disp_z;
  double p_disp_mag;
  double p_ke_rot;
  double p_ke_tra;
  int cell_x;
  int cell_y;
  int cell_z;
  long ts;
  std::string cellstr;
};

struct decomp_table {
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &cellstr;
  }
  std::string cellstr;
};

enum ContactTXTColumns {
  p1_x,
  p1_y,
  p1_z,
  p2_x,
  p2_y,
  p2_z,
  p1_vx,
  p1_vy,
  p1_vz,
  p2_vx,
  p2_vy,
  p2_vz,
  p1_id,
  p2_id,
  is_periodic,
  c_force_x,
  c_force_y,
  c_force_z,
  cn_force_x,
  cn_force_y,
  cn_force_z,
  ct_force_x,
  ct_force_y,
  ct_force_z,
  c_torque_x,
  c_torque_y,
  c_torque_z,
  disp_x,
  disp_y,
  disp_z,
  contact_area,
  contact_overlap,
  sliding_contact,
};

enum ParticleTXTColumns {
  id,
  type,
  x,
  y,
  z,
  rad,
  vx,
  vy,
  vz,
  fx,
  fy,
  fz,
  omegax,
  omegay,
  omegaz,
  coord,
  displace_x,
  displace_y,
  displace_z,
  disp_mag,
  ke_rot,
  ke_tra
};

enum dumpfile_t { particle, contact };

constexpr double Ystar = (6.5e11) / (2 * (1 - (.25 * .25)));
constexpr double Gstar = (6.5e11) / (4 * (2 - .25) * (1 + .25));

struct ts_column {
  long ts;
};

struct radius {
  int id;
  double rad;
};

using aggregate_map_per_particle =
    std::map<long, std::map<std::string, double>>;
using aggregate_map_t = std::map<std::string, std::map<std::string, double>>;

using cell_properties_map =
    std::map<std::string, std::map<long, std::map<std::string, double>>>;

struct agg_result_per_particle_t {
  std::string cell;
  aggregate_map_per_particle agg;
};

struct p_velocity_t {
  long ts;
  std::map<int, std::map<std::string, double>> vel;
};

struct aggr_result_t {
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &agg;
    ar &ts;
  }
  aggregate_map_t agg;
  long ts;
};

#endif
