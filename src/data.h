#ifndef DATA_H
#define DATA_H
#include <string>
#include <map>

typedef struct {
    std::string InputPath;
    std::string OutputPath;
} Config;


typedef struct{
    long ts;
    std::map<int,double> b_centrality;
    double mean;
} Result;
#endif


