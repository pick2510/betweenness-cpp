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

typedef struct {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &InputPath;
    ar &OutputPath;
    ar &sep;
  }

public:
  std::string InputPath;
  std::string OutputPath;
  char sep[2] = " ";
} Config;

typedef struct {
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

} Result;


constexpr char ts_particle_path[] = "particles";
constexpr char ts_centrality_path[] = "centrality";

#endif
