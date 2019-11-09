#include "compressor.h"
#include "partition.h"
#include "entropy_coding.h"
#include "elias_coding.h"

#include <sdsl/vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>

#include <iostream>

string compressor::from_binary(unsigned long long factor, const string &charset) {
  string res = "";
  unsigned long long base = charset.length()+1;
  while(factor) {
    unsigned long long div = factor % base;
    res = charset[div-1] + res;
    factor /= base;
  }

  return res;
}

unsigned long long compressor::to_binary(map<char, unsigned char> &alphabet_map, const string &s) {
  unsigned long long base = alphabet_map.size() + 1;
  unsigned long long result = 0;
  unsigned long long mul = 1;

  for(int i=int(s.length())-1; i>=0; --i) {
    unsigned long long letter = alphabet_map[s[i]] + 1;
    result += mul*letter;
    mul *= base;
  }
  return result;
}

std::tuple< vector<unsigned long long>, vector<unsigned int>, string>
  compressor::replace_partition_succinct(const vector<string> &partition) {

  //return values
  vector<unsigned long long> factors;
  vector<unsigned int> factor_orders;
  string charset = "";

  set<char> alphabet;
  for (unsigned int i=0; i<partition.size(); ++i)
    for(unsigned int j=0; j<partition[i].length(); ++j)
      alphabet.insert(partition[i][j]);

  assert(alphabet.size() < 256);


  map<char, unsigned char> old_new_char;
  unsigned char current = 0;
  for (auto it = alphabet.begin(); it != alphabet.end(); ++it) {
    old_new_char[*it] = current++;
    charset += *it;
  }

  set<string> dictionary;

  for (unsigned int i=0; i<partition.size(); ++i)
    dictionary.insert(partition[i]);

  vector< pair<unsigned long long, string> >  binary_to_factors;
  for (auto it = dictionary.begin(); it != dictionary.end(); ++it) {
    unsigned long long bin = to_binary(old_new_char, *it);
    binary_to_factors.push_back(make_pair(bin, *it));
  }
  std::sort(binary_to_factors.begin(), binary_to_factors.end());

  gp_hash_table<string, unsigned int> factors_to_order;

  for (unsigned int i=0; i<binary_to_factors.size(); ++i) {
    factors_to_order[binary_to_factors[i].second] = i;
    factors.push_back(binary_to_factors[i].first);
  }
  for (unsigned int i=0; i<partition.size(); ++i){
    factor_orders.push_back(factors_to_order[partition[i]]);
  }

  // for (uint i =0 ; i< factor_orders.size(); ++i )
  //   std::cout<<from_binary(factors[factor_orders[i]], charset );

  return std::make_tuple(factors, factor_orders, charset);
}

//Format for H0:
//size of dictionary(uint64)
// alphabet size (uchar), charset
//#of bits representing:
//  concatenated:
//    delta-compressed sequence of codes
//    delta-comperssed sequence of letters of replaced texts, corresponding to codes
//    delta compressed sequence of factors, corresponding to letters
//bitstring

vector<bool> compressor::compress_zeroth(int order, bool method, const string &s, bps_data* bps) {
  assert(order <= 8); //todo add support for longer factors

  vector<bool> result;
  auto part = method ?  quasi_optimal_partiton::do_partition(s, order, 1, -1)
                     :  quasi_optimal_partiton::best_shift_partition(0, s, order);
  auto replaced = replace_partition_succinct(part);
  auto codes = entropy_coding::get_huffman_codes_succinct(std::get<1>(replaced));

  unsigned long long dict_size = std::get<0>(codes).size();

  sdsl::int_vector<> codes_intv(dict_size, 0, 64);
  sdsl::int_vector<> letters_intv(dict_size, 0, 64);
  sdsl::int_vector<> factors_intv(dict_size, 0, 64);

  for (uint i=0; i<dict_size; ++i) {
    if (i == 0) {
      codes_intv[i] = std::get<0>(codes)[i];
      letters_intv[i] = std::get<1>(codes)[i];
      factors_intv[i] = std::get<0>(replaced)[i];
    } else {
      codes_intv[i] = std::get<0>(codes)[i] - std::get<0>(codes)[i-1];

      unsigned long long diff_abs =
        (std::get<1>(codes)[i] > std::get<1>(codes)[i-1])
        ? ( (std::get<1>(codes)[i] - std::get<1>(codes)[i-1]) << 1)
        : (( (std::get<1>(codes)[i-1] - std::get<1>(codes)[i]) << 1 )  + 1);
      letters_intv[i] = diff_abs;

      factors_intv[i] = std::get<0>(replaced)[i] - std::get<0>(replaced)[i-1];
    }
  }

  sdsl::int_vector<> concat_sequences(dict_size*3, 0, 64);
  for (uint i=0; i<dict_size; ++i) {
    concat_sequences[i] = codes_intv[i];
    concat_sequences[dict_size + i] = letters_intv[i];
    concat_sequences[dict_size*2 + i] = factors_intv[i];
  }

  // std::cout<<"Dict size: "<<dict_size<<"\n";

  //testing purposes
  // sdsl::int_vector<> codes_intv_encoded;
  // sdsl::int_vector<> letters_intv_encoded;
  // sdsl::int_vector<> factors_intv_encoded;
  //
  // sdsl::coder::elias_delta::encode(codes_intv, codes_intv_encoded);
  // sdsl::coder::elias_delta::encode(letters_intv, letters_intv_encoded);
  // sdsl::coder::elias_delta::encode(factors_intv, factors_intv_encoded);

  // std::cout<<"sizes: "<<sdsl::size_in_bytes(codes_intv_encoded)*8 << " "
  //          <<sdsl::size_in_bytes(letters_intv_encoded)*8 << " "
  //          <<sdsl::size_in_bytes(factors_intv_encoded) *8 <<" \n";

  sdsl::int_vector<> encoded_sequences;
  sdsl::coder::elias_delta::encode(concat_sequences, encoded_sequences);
  // std::cout<<"Total size in bits: "<<sdsl::size_in_bytes(encoded_sequences) *8 <<" \n";


  vector<bool> huffman_codes;
  for (uint i=0; i<(std::get<1>(replaced)).size(); ++i) {
    unsigned int letter = (std::get<1>(replaced))[i];
    string code = (std::get<2>(codes))[ letter ];
    for (uint j=0; j<code.length(); ++j)
      huffman_codes.push_back(code[j]-'0');
  }

  // std::cout<<"Entropy string size: "<<huffman_codes.size()<<"\n";



  //write to result

  vector<bool> empty64(64, 0);
  vector<bool> empty8(8, 0);
  *(empty64.begin()._M_p) = dict_size;
  result.insert(result.end(), empty64.begin(), empty64.end());

  string charset = std::get<2>(replaced);
  unsigned char alphabet_size = charset.length();
  *(empty8.begin()._M_p) = alphabet_size;

  result.insert(result.end(), empty8.begin(), empty8.end());
  for (uint i=0; i<alphabet_size; ++i) {
    *(empty8.begin()._M_p) = charset[i];
    result.insert(result.end(), empty8.begin(), empty8.end());
  }

  unsigned long long bits_to_encode_sequence = encoded_sequences.bit_size();
  *(empty64.begin()._M_p) = bits_to_encode_sequence;
  result.insert(result.end(), empty64.begin(), empty64.end());

  unsigned char* ptr_encoded_seq = (unsigned char*)(encoded_sequences.data());
  for (unsigned long long i=0; i<bits_to_encode_sequence;) {
    unsigned char buffer = *ptr_encoded_seq;
    *(empty8.begin()._M_p) = buffer;
    result.insert(result.end(), empty8.begin(), empty8.end());
    ptr_encoded_seq++;
    i += 8;
  }

  //bitstring
  unsigned long long huffman_bit_size = huffman_codes.size();
  *(empty64.begin()._M_p) = huffman_bit_size;
  result.insert(result.end(), empty64.begin(), empty64.end());

  result.insert(result.end(), huffman_codes.begin(), huffman_codes.end());

  if (bps != NULL) {
    bps->bps = double(result.size())/s.size();
    bps->dict_bps = double(result.size() - huffman_codes.size() )/s.size();
    bps->bitstring_bps = double(huffman_codes.size())/s.size();
  }

  return result;
}

string compressor::decompress_zeroth(const vector<bool> &v) {

  string s = "";
  unsigned char *ptr_read = (unsigned char*)(v.begin()._M_p);

  unsigned long long dict_size = (*(unsigned long long*)ptr_read);
  ptr_read += 8;
  unsigned char alphabet_size = (*(unsigned char*)ptr_read);
  ptr_read += 1;

  string charset = "";
  for (uint i=0; i < alphabet_size; ++i) {
    char symbol = (*(unsigned char*)ptr_read);
    charset += symbol;
    ptr_read += 1;
  }

  unsigned long long bits_to_decode_sequence = (*(unsigned long long*)ptr_read);
  ptr_read += 8;

  sdsl::int_vector<> encoded_sequences;
  encoded_sequences.bit_resize(bits_to_decode_sequence);
  encoded_sequences.width(64);

  unsigned char* ptr_encoded_seq = (unsigned char*)(encoded_sequences.data());
  for (unsigned long long i=0; i<bits_to_decode_sequence;) {
    *ptr_encoded_seq = *ptr_read;
    ptr_encoded_seq++;
    ptr_read++;
    i += 8;
  }

  sdsl::int_vector<> concat_sequences;
  sdsl::coder::elias_delta::decode(encoded_sequences, concat_sequences);

  vector<unsigned long long> codes;
  vector<unsigned long long> letters;
  vector<unsigned long long> factors;

  codes.push_back(concat_sequences[0]);
  letters.push_back(concat_sequences[dict_size]);
  factors.push_back(concat_sequences[dict_size*2]);

  for (unsigned int i=1; i<dict_size; ++i) {
    codes.push_back(codes[i-1] + concat_sequences[i]);

    unsigned long long letter = concat_sequences[dict_size + i] % 2 ?
        (letters[i-1] - (concat_sequences[dict_size + i] >> 1))
      : (letters[i-1] + (concat_sequences[dict_size + i] >> 1));
    letters.push_back(letter);

    factors.push_back(factors[i-1] + concat_sequences[dict_size*2 + i]);
  }

  // std::cout<<"Sequences read"<<std::endl;

  map<unsigned long long, unsigned long long> code_letter_map;
  for (unsigned int i=0; i<dict_size; ++i) {
    code_letter_map[codes[i]] = letters[i];
    // std::cout<<codes[i]<<" "<<entropy_coding::binary_to_code(codes[i])<<" "<<letters[i]<<"\n";
  }


  unsigned long long huffman_bit_size = (*(unsigned long long*)ptr_read);
  ptr_read += 8;

  vector<unsigned int> decoded_orders;
  unsigned char buffer = (*(unsigned char*)ptr_read);
  ptr_read++;
  unsigned char bits_in_buffer = 8;
  unsigned long long current_code = 1;
  // unsigned int current_code_length = 0;

  //TODO: add support for empty codes
  for (unsigned long long i=0; i<huffman_bit_size; ++i) {

    // if (i % 100000 == 0)
    //   std::cout<<"done "<< i << " out of "<<huffman_bit_size<<std::endl;

    if (!bits_in_buffer) {
      bits_in_buffer = 8;
      buffer = (*(unsigned char*)ptr_read);
      ptr_read++;
    }
    bool read_bit = buffer & (1<<(8-bits_in_buffer));
    current_code = ((current_code<<1) + read_bit);
    // current_code_length++;

    auto it = code_letter_map.find(current_code);
    if (it != code_letter_map.end()) {
      current_code = 1;
      decoded_orders.push_back(it->second);
    }

    bits_in_buffer--;
  }

  for (uint i =0 ; i< decoded_orders.size(); ++i )
    s += from_binary(factors[decoded_orders[i]], charset );

  return s;
}


void compressor::write_ull(vector<bool> &result, unsigned long long write_num) {
  vector<bool> empty64(64, 0);
  *(empty64.begin()._M_p) = write_num;
  result.insert(result.end(), empty64.begin(), empty64.end());
}

void compressor::read_ull(unsigned char* &ptr_read, unsigned long long &read_num) {
  read_num = (*(unsigned long long*)ptr_read);
  ptr_read += 8;
}

void compressor::write_byte(vector<bool> &result, unsigned char write_num) {
  vector<bool> empty8(8, 0);
  *(empty8.begin()._M_p) = write_num;
  result.insert(result.end(), empty8.begin(), empty8.end());
}

void compressor::read_byte(unsigned char* &ptr_read, unsigned char &read_num) {
  read_num = (*(unsigned char*)ptr_read);
  ptr_read++;
}


void compressor::write_elias_sequence(vector<bool> &result, const vector<unsigned long long> &write_vec) {

  sdsl::int_vector<> int_vec(write_vec.size(), 0, 64);
  for (unsigned long i=0; i<write_vec.size(); ++i)
    int_vec[i] = write_vec[i];

  sdsl::int_vector<> encoded_vec;
  sdsl::coder::elias_delta::encode(int_vec, encoded_vec);

  vector<bool> empty64(64, 0);
  vector<bool> empty8(8, 0);

  unsigned long long bits_to_encode = encoded_vec.bit_size();
  *(empty64.begin()._M_p) = bits_to_encode;
  result.insert(result.end(), empty64.begin(), empty64.end());

  unsigned char* ptr_encoded_seq = (unsigned char*)(encoded_vec.data());
  for (unsigned long long i=0; i<bits_to_encode;) {
    unsigned char buffer = *ptr_encoded_seq;
    *(empty8.begin()._M_p) = buffer;
    result.insert(result.end(), empty8.begin(), empty8.end());
    ptr_encoded_seq++;
    i += 8;
  }
}

void compressor::read_elias_sequence(unsigned char* &ptr_read, vector<unsigned long long> &read_vec) {

  unsigned long long bits_to_decode_sequence = (*(unsigned long long*)ptr_read);
  ptr_read += 8;

  sdsl::int_vector<> encoded_sequences;
  encoded_sequences.bit_resize(bits_to_decode_sequence);
  encoded_sequences.width(64);

  unsigned char* ptr_encoded_seq = (unsigned char*)(encoded_sequences.data());
  for (unsigned long long i=0; i<bits_to_decode_sequence;) {
    *ptr_encoded_seq = *ptr_read;
    ptr_encoded_seq++;
    ptr_read++;
    i += 8;
  }

  read_vec.clear();

  sdsl::int_vector<> decoded_sequences;
  sdsl::coder::elias_delta::decode(encoded_sequences, decoded_sequences);

  for(unsigned long i=0; i < decoded_sequences.size(); ++i)
    read_vec.push_back(decoded_sequences[i]);
}


//Format for H1:
//#of dictionaries
//elias compressed sizes of dictionaries
// alphabet size (uchar), charset
//#of bits representing:
//  concatenated:
//    delta-compressed sequence of codes
//    delta-comperssed sequence of letters of replaced texts, corresponding to codes
// #of different factors
//  delta compressed sequence of factors, corresponding to letters
//bitstring

vector<bool> compressor::compress_first(int order, bool method, const string &s, bps_data* bps) {
  assert(order <= 4); //todo add support for longer factors

  vector<bool> result;
  auto part = method ?  quasi_optimal_partiton::do_partition_first(s, order, 1, -1) //limit 10^5
                     :  quasi_optimal_partiton::best_shift_partition(1, s, order);
  auto replaced = replace_partition_succinct(part);
  auto codes = entropy_coding::get_huffman_codes_succinct_first_order(std::get<1>(replaced));


  unsigned long long different_factors = codes.size();
  write_ull(result, different_factors);

  vector<unsigned long long> dict_sizes;
  for (unsigned long i=0; i<codes.size(); ++i)
    dict_sizes.push_back((unsigned long long)(std::get<2>(codes[i]).size()));

  write_elias_sequence(result, dict_sizes);

  string charset = std::get<2>(replaced);
  unsigned char alphabet_size = charset.length();
  write_byte(result, alphabet_size);

  for (uint i=0; i<alphabet_size; ++i)
    write_byte(result, charset[i]);

  vector<unsigned long long> codes_delta;
  vector<unsigned long long> letters_delta;

  unsigned long long prev_letter = 0;
  for (unsigned long i=0; i<codes.size(); ++i) {
    //encode dictionary
    for (unsigned long j=0; j<dict_sizes[i]; ++j){
      unsigned long long code =  std::get<0>(codes[i])[j];
      if (j == 0)
        codes_delta.push_back(code);
      else
        codes_delta.push_back(code - std::get<0>(codes[i])[j-1]);

      unsigned long long current_letter = std::get<1>(codes[i])[j];
      if (i == 0 && j == 0 ){
        letters_delta.push_back(current_letter);
      } else {
        unsigned long long diff_abs =
          (current_letter > prev_letter)
          ? ( (current_letter - prev_letter) << 1)
          : (( (prev_letter - current_letter) << 1 )  + 1);
        letters_delta.push_back(diff_abs);
      }
      prev_letter = current_letter;
    }
  }

  write_elias_sequence(result, codes_delta);
  write_elias_sequence(result, letters_delta);

  vector<unsigned long long> factors_delta;
  factors_delta.push_back(std::get<0>(replaced)[0]);
  for (unsigned long i=1; i<std::get<0>(replaced).size(); ++i) {
    factors_delta.push_back(std::get<0>(replaced)[i] - std::get<0>(replaced)[i-1]);
  }
  write_elias_sequence(result, factors_delta);

  unsigned long long symbols_to_write = (std::get<1>(replaced)).size();
  write_ull(result, symbols_to_write);

  unsigned long long first_letter = (std::get<1>(replaced))[0];
  write_ull(result, first_letter);

  vector<bool> huffman_codes;
  unsigned int context_letter = first_letter;
  for (uint i=1; i<(std::get<1>(replaced)).size(); ++i) {
    unsigned int letter = (std::get<1>(replaced))[i];
    string code = (std::get<2>(codes[context_letter]))[ letter ];
    for (uint j=0; j<code.length(); ++j)
      huffman_codes.push_back(code[j]-'0');
    context_letter = letter;
  }

  // unsigned long long huffman_bit_size = huffman_codes.size();
  // write_ull(result, huffman_bit_size);

  result.insert(result.end(), huffman_codes.begin(), huffman_codes.end());

  // std::cout<< "total: "<<result.size()
  //          << " dict: "<<result.size() -  huffman_codes.size()
  //          <<" huffman: "<<huffman_codes.size()<<"\n";

  // std::cout<<first_letter<<" fl\n";

  if (bps != NULL) {
    bps->bps = double(result.size())/s.size();
    bps->dict_bps = double(result.size() - huffman_codes.size() )/s.size();
    bps->bitstring_bps = double(huffman_codes.size())/s.size();
  }


  return result;
}


string compressor::decompress_first(const vector<bool> &v) {
  string s;

  unsigned char *ptr_read = (unsigned char*)(v.begin()._M_p);

  unsigned long long different_factors;
  read_ull(ptr_read, different_factors);

  vector<unsigned long long> dict_sizes;
  read_elias_sequence(ptr_read, dict_sizes);

  unsigned char alphabet_size;
  read_byte(ptr_read, alphabet_size);
  string charset = "";
  for (uint i=0; i<alphabet_size; ++i) {
    unsigned char letter;
    read_byte(ptr_read, letter);
    charset += letter;
  }

  vector<unsigned long long> codes_delta;
  vector<unsigned long long> letters_delta;

  read_elias_sequence(ptr_read, codes_delta);
  read_elias_sequence(ptr_read, letters_delta);

  vector<vector<unsigned long long>> codes_in_dict(different_factors);
  vector<vector<unsigned long long>> letters_in_dict(different_factors);


  unsigned long codes_delta_index = 0;
  unsigned long long current_letter = 0;
  for (unsigned long i=0; i<different_factors; ++i) {
    //encode dictionary
    unsigned long long current_code = 0;
    for (unsigned long j=0; j<dict_sizes[i]; ++j){
      if (j == 0)
        codes_in_dict[i].push_back(codes_delta[codes_delta_index]);
      else
        codes_in_dict[i].push_back(current_code + codes_delta[codes_delta_index]);
      current_code += codes_delta[codes_delta_index];

      if (i == 0  && j == 0)
        current_letter = letters_delta[codes_delta_index];
      else {
        current_letter = letters_delta[codes_delta_index] % 2 ?
            (current_letter - (letters_delta[codes_delta_index] >> 1))
          : (current_letter + (letters_delta[codes_delta_index] >> 1));
      }
      letters_in_dict[i].push_back(current_letter);

      codes_delta_index++;
    }
  }

  vector<map<unsigned long long, unsigned long long>> code_letter_map(different_factors);
  for (unsigned long i=0; i<different_factors; ++i) {
    for (unsigned long j=0; j<dict_sizes[i]; ++j)
        code_letter_map[i][codes_in_dict[i][j]] = letters_in_dict[i][j];
  }

  vector<unsigned long long> factors_delta;
  read_elias_sequence(ptr_read, factors_delta);
  vector<unsigned long long> factors;
  factors.push_back(factors_delta[0]);
  for (unsigned long i = 1; i<factors_delta.size(); ++i) {
    factors.push_back(factors.back() + factors_delta[i]);
  }


  unsigned long long symbols_to_read;
  read_ull(ptr_read, symbols_to_read);


  unsigned long long first_letter;
  read_ull(ptr_read, first_letter);

  // unsigned long long huffman_bit_size;
  // read_ull(ptr_read, huffman_bit_size);

  vector<unsigned int> decoded_orders;
  decoded_orders.push_back(first_letter);



  unsigned char buffer = (*(unsigned char*)ptr_read);
  ptr_read++;
  unsigned char bits_in_buffer = 8;

  unsigned long context = first_letter;
  unsigned long long current_code = 1;

  for (unsigned long long read_symbols = 1; read_symbols < symbols_to_read; ++read_symbols) {
      auto it = code_letter_map[context].find(current_code);
      while (it == code_letter_map[context].end()) {
        if (!bits_in_buffer) {
          bits_in_buffer = 8;
          buffer = (*(unsigned char*)ptr_read);
          ptr_read++;
        }
        bool read_bit = buffer & (1<<(8-bits_in_buffer));
        current_code = ((current_code<<1) + read_bit);
        bits_in_buffer--;
        it = code_letter_map[context].find(current_code);
      }

      current_code = 1;
      decoded_orders.push_back(it->second);
      context = it->second;
  }

  for (uint i =0 ; i< decoded_orders.size(); ++i )
    s += from_binary(factors[decoded_orders[i]], charset );


  return s;
}
