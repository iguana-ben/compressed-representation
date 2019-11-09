#ifndef ENTROPY_UTILS_H
#define ENTROPY_UTILS_H

#include <vector>
#include <map>
#include <string>

using std::vector;
using std::map;
using std::string;


double bit_size(int order, const vector<int> &str);
double bit_size(int order, const string &str);

double h0_vec(const vector<int> &str);
double h1_vec(const vector<int> &str);


int count_fo_dictionary_sizes(const vector<int> &str);
unsigned long long count_fo_contexts(const vector<int> &str);

#endif
