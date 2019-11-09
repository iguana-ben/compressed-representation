#ifndef COMPRESSED_ARRAY_H
#define COMPRESSED_ARRAY_H

#include <vector>
// #include <sdsl/bit_vector.hpp>
#include <sdsl/vectors.hpp>

#include "compressed_vector_rrr.h"
#include "compressor.h"

using std::vector;

#define MAX_ALPHABET_SIZE 256
#define MAX_CODE_LENGTH 58 //64 - log(64)


// class simple_prefix_sum {
// private:
//   unsigned long long top_level_div;
//   sdsl::enc_vector<> encoded_vec_top_level;
//   sdsl::enc_vector<> encoded_vec;
// public:
//   simple_prefix_sum(const vector<unsigned int> &v);
//   unsigned long long get_size_in_bits();
//   unsigned long long find(unsigned long long i);
//   unsigned long long get(unsigned long long i);
// }

class compressed_array_zeroth {
  bool method;
  compressed_vector_rrr* compressed_codes_rrr;

  unsigned long long max_factor_length_;
  unsigned long long code_len_mn;
  compressed_vector_rrr* compressed_borders_rrr;

  vector<bool> huffman_codes;

  string charset;

  sdsl::enc_vector<> encoded_factors;
  sdsl::enc_vector<> encoded_factors_in_dict;

  unsigned long long mod_lookup;
  sdsl::enc_vector<> encoded_lengths_lookup;
  sdsl::int_vector<> lengths_offsets_intv;
  vector<unsigned long long> lengths_precompute;
  // vector<unsigned long long> lengths_lookup;
  // vector<unsigned long long> lengths_offsets;


  int first_factor_len;

  vector<unsigned long long> precompute_powers;

  unsigned long long read_factor(const unsigned long long &code_start, const unsigned long long &code_end);
  inline unsigned char get_length_from_factor(const unsigned long long &factor);
  // inline unsigned char get_length_from_factor(const unsigned long long &factor);
  inline void fast_decode_factor(char *buffer, const unsigned long long &factor);
  inline void fast_decode_factor_w_len(char *buffer, const unsigned long long &factor, unsigned char &factor_len);
public:
  compressed_array_zeroth(string &v, int max_factor_length, bool length_penalty, bool method_);
  ~compressed_array_zeroth() {
    delete compressed_codes_rrr;
    delete compressed_borders_rrr;
  }
  char get(int i);
  void read_block(char* out_buffer, unsigned long start, unsigned long num );

  unsigned long long get_size_in_bits() {
    unsigned long long sm = compressed_codes_rrr->get_size_in_bytes()*8 +
                            compressed_borders_rrr->get_size_in_bytes()*8 +
                            64*5 + 8*charset.length() +
                            huffman_codes.size() +
                            sdsl::size_in_bytes(encoded_factors)*8 +
                            sdsl::size_in_bytes(encoded_factors_in_dict)*8;
    if (method)
      sm += sdsl::size_in_bytes(encoded_lengths_lookup)*8 +
            sdsl::size_in_bytes(lengths_offsets_intv)*8 + 8;
            //+ lengths_precompute.size()*64;

    return sm;
  }
};

#endif
