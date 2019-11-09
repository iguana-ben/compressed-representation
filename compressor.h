#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <vector>
#include <string>
#include <set>
#include <map>
#include <ext/pb_ds/assoc_container.hpp>


using std::vector;
using std::string;
using std::set;
using std::map;
using namespace __gnu_pbds;


class compressor {
private:
  static unsigned long long to_binary(map<char, unsigned char> &alphabet_map, const string &partition);

  static void write_ull(vector<bool> &result, unsigned long long write_num);
  static void read_ull(unsigned char* &ptr_read, unsigned long long &read_num);

  static void write_byte(vector<bool> &result, unsigned char write_num);
  static void read_byte(unsigned char* &ptr_read, unsigned char &read_num);

  static void write_elias_sequence(vector<bool> &result, const vector<unsigned long long> &write_vec);
  static void read_elias_sequence(unsigned char* &ptr_read, vector<unsigned long long> &read_vec);

public:
  struct bps_data {
    double bps;
    double dict_bps;
    double bitstring_bps;
  };

  static string from_binary(unsigned long long factor, const string &charset);
  static std::tuple< vector<unsigned long long>, vector<unsigned int>, string >
    replace_partition_succinct(const vector<string> &s);

  static vector<bool> compress_zeroth(int order, bool method, const string &s, bps_data* bps);
  static string decompress_zeroth(const vector<bool> &v);

  static vector<bool> compress_first(int order, bool method, const string &s, bps_data* bps);
  static string decompress_first(const vector<bool> &v);
};

#endif
