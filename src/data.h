#ifndef DATA_H
#define DATA_H
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>


typedef struct {
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
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


typedef struct{
private: 
   friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar &ts;
        ar &b_centrality;
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
    std::map<int,double> b_centrality;
    double mean;
    double var;
    double std;
    double skew;
    double kur;
    double q_090;
    double q_099;
  
} Result;


#endif


