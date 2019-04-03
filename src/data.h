#ifndef DATA_H
#define DATA_H
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>



typedef struct {
    std::string InputPath;
    std::string OutputPath;
    char sep[2] = " ";
} Config;


typedef struct{
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


