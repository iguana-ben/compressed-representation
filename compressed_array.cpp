#include "compressed_array.h"

#include <sdsl/vectors.hpp>
#include <ctime>
#include <bitset>


#include "entropy_coding.h"
#include "partition.h"

// simple_prefix_sum::simple_prefix_sum(const vector<unsigned int> &v) {
//   top_level_div = 64;
//   sdsl::int_vector<> int_vec(factors.size(), 0, 64);
//   for (unsigned long i=0; i<factors.size(); ++i)
//     int_vec[i] = factors[i];
// }

compressed_array_zeroth::compressed_array_zeroth(string &v, int max_factor_length,
  bool length_penalty, bool method_) {
    assert(max_factor_length <= 8); //todo: add support for longer factors

    max_factor_length_ = max_factor_length;
    method = method_;
    auto part = method ?
      quasi_optimal_partiton::do_partition(v, max_factor_length, length_penalty, -1):
      quasi_optimal_partiton::best_shift_partition(0, v, max_factor_length);

    vector<unsigned long long> factors;
    vector<unsigned int> replaced_sequence;
    std::tie(factors, replaced_sequence, charset) = compressor::replace_partition_succinct(part);

    vector<unsigned long long> codes;
    vector<unsigned int> factors_in_dict;
    map<unsigned int, string> codes_map;
    std::tie(codes, factors_in_dict, codes_map) = entropy_coding::get_huffman_codes_succinct(replaced_sequence);

    // for (unsigned long i=0; i<factors_in_dict.size(); ++i) {
    //   std::cout<<codes[i]<<" "<<factors[i]<<" "<<factors_in_dict[i]<<"\n";
    // }

    sdsl::int_vector<> codes_intv(codes.size(), 0, 64);
    sdsl::int_vector<> factors_intv(factors.size(), 0, 64);
    sdsl::int_vector<> factors_in_dict_intv(factors_in_dict.size(), 0, 64);

    vector<unsigned long long> codes_delta(codes.size());
    codes_delta[0] = codes[0];
    for (unsigned long i=0; i<codes.size(); ++i) {
      codes_intv[i] = codes[i];
      if (i)
        codes_delta[i] = codes[i]-codes[i-1];
    }

    for (unsigned long i=0; i<factors.size(); ++i)
      factors_intv[i] = factors[i];

    for (unsigned long i=0; i<factors_in_dict.size(); ++i)
      factors_in_dict_intv[i] = factors_in_dict[i];

    compressed_codes_rrr = new compressed_vector_rrr(codes_delta);
    // std::cout<<" compressed with rrr: "<<compressed_codes_rrr->get_size_in_bytes() * 8 <<"\n";

    // auto compressed_factors_sd = new coxmpressed_vector_sd(factors);
    // std::cout<<" compressed with sd: "<<compressed_factors_sd->get_size_in_bytes() * 8 <<"\n";



    sdsl::enc_vector<> encoded_codes(codes_intv);
    encoded_factors = sdsl::enc_vector<>(factors_intv);
    encoded_factors_in_dict = sdsl::enc_vector<>(factors_in_dict_intv);

    // std::cout<<sdsl::size_in_bytes(encoded_codes)*8 <<" "
    //          <<sdsl::size_in_bytes(encoded_factors)*8 <<" "
    //          <<sdsl::size_in_bytes(encoded_factors_in_dict)*8 <<"\n";


     huffman_codes = vector<bool>();

     vector<unsigned long long> code_borders;

     vector<unsigned long long> factor_lengths;
     unsigned long long dv = max_factor_length;
     unsigned long long sm = 0;

     code_len_mn = 64;

     vector<unsigned long long> lengths_lookup;
     vector<unsigned long long> lengths_offsets;
     // lengths_lookup.clear();
     // lengths_offsets.clear();

     lengths_lookup.push_back(0);
     lengths_offsets.push_back(0);
     unsigned long long sum_lengths = 0;
     mod_lookup = max_factor_length*8;

     for (unsigned long i=0; i<replaced_sequence.size(); ++i) {
       unsigned int letter = replaced_sequence[i];
       string code = codes_map[letter];
       for (uint j=0; j<code.length(); ++j)
         huffman_codes.push_back(code[j]-'0');

       code_borders.push_back(code.length());

       // if( i < 1)
       //  std::cout<<code << " "<< letter<<" code\n";

       // std::cout<<code.length() <<"\n";

       if(code.length() < code_len_mn)
        code_len_mn = code.length();


       string factor = compressor::from_binary(factors[letter], charset);
       if (i == 0 ) first_factor_len = factor.length();
       // std::cout<<factor<<"\n";
       // if (factor.length() != (unsigned int)max_factor_length)
       sm += factor.length();
       factor_lengths.push_back(sm);
       if (i % dv == 0) {
        sm = 0;
       }
       // 63 64 65 66 67
       // x y
       // 67
       if (sum_lengths/mod_lookup != (sum_lengths + factor.length())/mod_lookup) {
         lengths_lookup.push_back(i);
         unsigned long long closest_mul = (sum_lengths + factor.length()) - (sum_lengths + factor.length()) %mod_lookup;
         // closest_mul *= mod_lookup;
         // sum_lengths - sum_lengths/mod_lookup
         lengths_offsets.push_back(closest_mul - sum_lengths);
         // std::cout<<i << " "<<(sum_lengths + factor.length()) % mod_lookup << " lo\n";
       }
       sum_lengths += factor.length();



       // std::cout<< factor << " "<< factor.length() << std::endl;

     }


     sdsl::int_vector<> lengths_lookup_intv(lengths_lookup.size(), 0, 64);
     lengths_offsets_intv = sdsl::int_vector<>(lengths_offsets.size(), 0, 64);

     vector<unsigned long long> lengths_lookup_differences(lengths_lookup.size());
     lengths_lookup_differences[0] = lengths_lookup[0];

     for(unsigned long i=1; i<lengths_lookup.size(); ++i)
        lengths_lookup_differences[i] = lengths_lookup[i] - lengths_lookup[i-1];

     for(unsigned long i=0; i<lengths_lookup.size(); ++i) {
      lengths_lookup_intv[i] = lengths_lookup[i];
      lengths_offsets_intv[i] = lengths_offsets[i];
    }

     encoded_lengths_lookup = sdsl::enc_vector<>(lengths_lookup_intv);
     sdsl::util::bit_compress(lengths_offsets_intv);


     auto compressed_lengths_probe_rrr = compressed_vector_rrr(lengths_lookup_differences);
     auto compressed_lengths_rrr = compressed_vector_rrr(factor_lengths);
     // std::cout<<"compressed lengths structure: " << sdsl::size_in_bytes(encoded_lengths_lookup)*8 <<"\n";
     // std::cout<<" compressed with rrr: "<<compressed_lengths_probe_rrr.get_size_in_bytes() * 8 <<"\n";
     // std::cout<<" offsets:: "<<sdsl::size_in_bytes(lengths_offsets_intv) * 8 <<"\n";
     // std::cout<<" compressed with rrr: "<<compressed_lengths_rrr.get_size_in_bytes() * 8 <<"\n";




     for (unsigned long i = 0;i < code_borders.size(); ++i )
      code_borders[i] -= code_len_mn;


     compressed_borders_rrr = new compressed_vector_rrr(code_borders);
     // std::cout<<" compressed with rrr borders: "<<compressed_borders_rrr->get_size_in_bytes() * 8 <<"\n";


     // std::cout<<replaced_sequence.size() << " \n";
     // std::cout<<huffman_codes.size() << " \n";

    // clock_t begin = clock();
    // for (unsigned int i=0; i< 100000; ++i) {
    //     std::cout<<" 0th code "<<compressed_borders_rrr->sum(replaced_sequence.size()-i-5) <<" "
    //          << compressed_borders_rrr-> sum(replaced_sequence.size()-i-4) <<"\n";
    //     // std::cout<<" 0th code "<<code_borders[replaced_sequence.size()-i-5] <<" "
    //              // << code_borders[replaced_sequence.size()-i-4] <<"\n";
    //  }
     // clock_t end = clock();
     // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
     // std::cout<<elapsed_secs << " \n";

     // std::cout<<"len compressedr rrr: " << compressed_lengths_rrr .get_size_in_bytes() * 8<< "\n";

     // unsigned long long total_bits = huffman_codes.size() + compressed_borders_rrr->get_size_in_bytes() * 8 +
     //                                                        compressed_codes_rrr->get_size_in_bytes()*8 +
     //                                                        sdsl::size_in_bytes(encoded_factors)*8 +
     //                                                        sdsl::size_in_bytes(encoded_factors_in_dict)*8 +
     //                                                        1*(compressed_lengths_rrr.get_size_in_bytes() * 8);

    // std::cout<<"total bits: "<<total_bits <<" bps: " << double(total_bits)/v.length() <<std::endl;

    // for(uint i=0; i<128; ++i )
    //   std::cout<<huffman_codes[i];
    // std::cout<<"\n";

    unsigned long long to_add = 64 + huffman_codes.size() % 64; //padding
    while(to_add--) huffman_codes.push_back(0);
    for (unsigned long i=0; i<huffman_codes.size()/64; ++i)
      for(unsigned long j=0; j<32; ++j) {
        std::swap(huffman_codes[i*64 +j], huffman_codes[i*64 + 64-j-1]);
      }

    // for(uint i=0; i<128; ++i )
    //   std::cout<<huffman_codes[i];
    // std::cout<<"\n";

    unsigned long long base = charset.length()+1 ;
    unsigned long long precompute_len_min = 1;
    // unsigned long long precompute_len_max = (base-1);
    // lengths_precompute = vector<unsigned long long>(64);
    lengths_precompute.clear();
    // lengths_precompute.push_back(0);
    for (uint i=0; i<max_factor_length_; ++i) {
      // std::cout<<precompute_len_min<<" "<<64-__builtin_clzll(precompute_len_min)<< "  pl\n";
      // std::cout<<precompute_len_max<<" "<<64-__builtin_clzll(precompute_len_max)<< "  pl\n\n";
      // int len_min = 64-__builtin_clzll(precompute_len_min);
      // int len_max = 64-__builtin_clzll(precompute_len_max);
      lengths_precompute.push_back(precompute_len_min);
      precompute_len_min *= base;
      // precompute_len_max = precompute_len_max*base + (base-1);
    }

    precompute_powers.clear();
    precompute_powers.push_back(1);
    for (uint i=0; i<max_factor_length_; ++i)
      precompute_powers.push_back(precompute_powers.back()*base);

}

unsigned long long compressed_array_zeroth::read_factor(
  const unsigned long long &code_start, const unsigned long long &code_end) {

  unsigned long long shift_ptr =  (code_start>>6);
  unsigned long long* data_ptr = (unsigned long long*)huffman_codes.begin()._M_p;
  unsigned long long first_chunk = *(data_ptr + (shift_ptr));

  unsigned long long factor = 0;
  if ((shift_ptr+1)*64 > code_end ) {
    unsigned long long mask = -1;
    mask >>= (code_start & 255);
    factor = ((first_chunk & mask)>>(63-code_end)) + (1<<(code_end-code_start+1));
  } else {
    unsigned long long second_chunk = *(data_ptr + (shift_ptr + 1));
    unsigned long long mask = -1;
    mask >>= (code_start & 255);
    factor = ((first_chunk & mask)<<((code_end & 255) +1)) + (second_chunk>>(63- (code_end & 255)));
    factor += (1<<(code_end-code_start+1));
  }

  unsigned long long index=compressed_codes_rrr->search(factor);
  return encoded_factors[encoded_factors_in_dict[index]];
};

unsigned char compressed_array_zeroth::get_length_from_factor(const unsigned long long &factor) {
  // unsigned char result = max_factor_length_;
  // while(lengths_precompute[result-1] > factor) result--;
  unsigned long long base = charset.length() +1;
  unsigned char result = 0;
  unsigned long long factor_temp = factor;
  while(factor_temp > 0){
    result ++;
    factor_temp /= base;
  }
  return result;
}

void compressed_array_zeroth::fast_decode_factor(char *buffer, const unsigned long long &factor) {
  unsigned char index = 0;
  unsigned long long factor_temp = factor;
  unsigned long long div = factor;
  unsigned long long base = charset.size() + 1;

  while ( factor_temp ) {
    div /= base;
    buffer[index++] = charset[factor_temp -  div*base - 1];
    factor_temp = div;
  }
  unsigned char i = 0;
  index--;
  while(i < index)
    std::swap(buffer[i++], buffer[index--]);
}

void compressed_array_zeroth::fast_decode_factor_w_len(char *buffer,
    const unsigned long long &factor,
    unsigned char &factor_len) {
  unsigned char index = 0;
  unsigned long long factor_temp = factor;
  unsigned long long div = factor;
  unsigned long long base = charset.size() + 1;

  while ( factor_temp ) {
    div /= base;
    buffer[index++] = charset[factor_temp -  div*base - 1];
    factor_temp = div;
  }
  factor_len = index;

  unsigned char i = 0;
  index--;
  while(i < index)
    std::swap(buffer[i++], buffer[index--]);
}

char compressed_array_zeroth::get(int i) {
  char buffer[8];
  memset(buffer, 0, sizeof(buffer));
  if (method) {
    unsigned long long dv_i =  i/mod_lookup;
    unsigned long long block_index = encoded_lengths_lookup[dv_i];
    unsigned long long offset = lengths_offsets_intv[dv_i];
    unsigned long long dist = i - (dv_i*mod_lookup - offset);
    unsigned long long factor;
    unsigned long long code_start =
      compressed_borders_rrr->sum(block_index)+ block_index*code_len_mn;
    while(true) {

      // unsigned long long code_end =
        // compressed_borders_rrr->sum(block_index+1)-1 + (block_index+1)*code_len_mn;
      unsigned long long bit_index = code_start - block_index*code_len_mn + block_index;
      unsigned long long bit_array =  compressed_borders_rrr->get_next_64(bit_index);
      unsigned long long code_end = code_start + __builtin_ctzll(bit_array) + code_len_mn-1;

      factor = read_factor(code_start, code_end);
      unsigned char factor_length = get_length_from_factor(factor);
      if (dist < factor_length)
        break;

      dist -= factor_length;
      block_index++;

      code_start = code_end + 1;
    }
    fast_decode_factor(buffer, factor);
    return buffer[dist];
  } else {
    unsigned long long block_index = i < first_factor_len ? 0
                  :((i-first_factor_len)/max_factor_length_ + 1);
    unsigned long long in_block_index = block_index == 0 ? i : i - (first_factor_len + (block_index-1)*max_factor_length_);
    unsigned long long code_start =
      compressed_borders_rrr->sum(block_index)+ block_index*code_len_mn;

    // unsigned long long code_end =
      // compressed_borders_rrr->sum(block_index+1)-1 + (block_index+1)*code_len_mn;
    unsigned long long bit_index = code_start - block_index*code_len_mn + block_index;
    unsigned long long bit_array =  compressed_borders_rrr->get_next_64(bit_index);
    unsigned long long code_end = code_start + __builtin_ctzll(bit_array) + code_len_mn-1;

    unsigned long long factor = read_factor(code_start, code_end);
    fast_decode_factor(buffer, factor);

    return buffer[in_block_index];
  }

  return 0;
}

void compressed_array_zeroth::read_block(char* out_buffer, unsigned long start, unsigned long num )  {
  char buffer[8];
  memset(buffer, 0, sizeof(buffer));

  unsigned long long factor = 0;

  unsigned long long code_start = 0;
  unsigned long long code_end = 0;

  unsigned long long block_index = 0;
  unsigned long long dist = 0;
  unsigned char factor_length = 0;

  if (method) {
    unsigned long long dv_i =  start/mod_lookup;
    block_index = encoded_lengths_lookup[dv_i];
    unsigned long long offset = lengths_offsets_intv[dv_i];
    dist = start - (dv_i*mod_lookup - offset);

    code_start = compressed_borders_rrr->sum(block_index)+ block_index*code_len_mn;

    while(true) {

      unsigned long long bit_index = code_start - block_index*code_len_mn + block_index;
      unsigned long long bit_array =  compressed_borders_rrr->get_next_64(bit_index);
      code_end = code_start + __builtin_ctzll(bit_array) + code_len_mn-1;
      // code_end = compressed_borders_rrr->sum(block_index+1)-1 + (block_index+1)*code_len_mn;
      factor = read_factor(code_start, code_end);
      factor_length = get_length_from_factor(factor);
      if (dist < factor_length)
        break;

      dist -= factor_length;
      block_index++;
      code_start = code_end + 1;
    }
  } else {
    block_index = start < (unsigned long)first_factor_len ? 0
                  :((start-first_factor_len)/max_factor_length_ + 1);

    code_start = compressed_borders_rrr->sum(block_index)+ block_index*code_len_mn;

    unsigned long long bit_index = code_start - block_index*code_len_mn + block_index;
    unsigned long long bit_array =  compressed_borders_rrr->get_next_64(bit_index);
    code_end = code_start + __builtin_ctzll(bit_array) + code_len_mn-1;
    // code_end = compressed_borders_rrr->sum(block_index+1)-1 + (block_index+1)*code_len_mn;

    dist = block_index == 0 ? start : start - (first_factor_len + (block_index-1)*max_factor_length_);
    factor = read_factor(code_start, code_end);
    factor_length = block_index == 0 ? first_factor_len : max_factor_length_;
  }
  fast_decode_factor(buffer, factor);

  unsigned long long read_cnt = 0;
  block_index++;
  // code_start = code_end + 1;
  while(read_cnt != num) {
    if (dist == factor_length) {
      code_start = code_end + 1;
      unsigned long long bit_index = code_start - block_index*code_len_mn + block_index;
      unsigned long long bit_array =  compressed_borders_rrr->get_next_64(bit_index);
      code_end = code_start + __builtin_ctzll(bit_array) + code_len_mn-1;
      // code_end = compressed_borders_rrr->sum(block_index+1)-1 + (block_index+1)*code_len_mn;
      block_index++;
      factor = read_factor(code_start, code_end);

      if (method) {
        fast_decode_factor_w_len(buffer, factor, factor_length);
      } else {
        fast_decode_factor(buffer, factor);
        factor_length = max_factor_length_;
      }

      dist = 0;
    }
    out_buffer[read_cnt++] = buffer[dist++];
  }
}
