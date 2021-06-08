//this component was written by Jeffrey P. Spence in 2016
// email: spence.jeffrey@berkeley.edu

#ifndef LDHELMET_CONVERT_TABLE_CONVERT_TABLE_COMPONENT_H_
#define LDHELMET_CONVERT_TABLE_CONVERT_TABLE_COMPONENT_H_

#include <string>
#include <unordered_map>
#include "common/read_confs.h"
int ConvertTableComponent(std::string const base_command, int argc, char **argv);
bool isSubconfig(Conf sub_config, Conf super_config);
std::vector<double> marginalProbs(Conf config, std::unordered_map<std::string, std::vector<double>> &like_map);
std::vector<double> vectorSum(std::vector<double> vec1, std::vector<double> vec2);
static std::vector<double> logfactorials;
#endif // LDHELMET_CONVERT_TABLE_CONVERT_TABLE_H_
