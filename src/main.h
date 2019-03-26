#ifndef MAIN_H
#define MAIN_H

#include <string>
typedef struct {
    std::string InputPath;
    std::string OutputPath;
} Config;

Config getCL(int &arg, char **argv);

#endif