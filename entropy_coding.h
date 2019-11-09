#ifndef ENTROPY_CODING_H
#define ENTROPY_CODING_H

#include <map>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <set>



using std::map;
using std::string;
using std::vector;
using std::pair;

map<unsigned int, string> gen_huffman_codes(const vector<unsigned int> &text);
// pair<vector<unsigned long long>, map<string, int>> map_to_numbers(const vector<string> &v);

class entropy_coding {
private:
  static unsigned long long code_to_binary(const string &code);
  // static vector<bool>

public:
  static string binary_to_code(unsigned long long num);
  static std::tuple<vector<unsigned long long>, vector<unsigned int>, map<unsigned int, string> >
    get_huffman_codes_succinct(const vector<unsigned int> &text);
  static vector<std::tuple<vector<unsigned long long>, vector<unsigned int>, map<unsigned int, string> > >
    get_huffman_codes_succinct_first_order(const vector<unsigned int> &text);
};
#endif
