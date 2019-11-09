#include "compressed_vector_rrr.h"

#include <bitset>

compressed_vector_rrr::compressed_vector_rrr(vector<unsigned long long> &v) {

  unsigned long long sm = 0;
  for (auto it = v.begin(); it != v.end(); ++it)
    sm += (*it);
  sm += 66; //padding for get_next_64
  sdsl::bit_vector bv(sm + v.size() - 1, 0);

  sm = 0;
  for (unsigned long long i=0; i<v.size(); ++i) {
    sm += v[i];
    bv[i+sm] = 1;
  }
  rv = sdsl::rrr_vector<>(bv);


  rv_rank_0 = sdsl::rrr_vector<>::rank_0_type(&rv);
  rv_rank_1 = sdsl::rrr_vector<>::rank_1_type(&rv);
  rv_select_0 = sdsl::rrr_vector<>::select_0_type(&rv);
  rv_select_1 = sdsl::rrr_vector<>::select_1_type(&rv);

}


unsigned long long compressed_vector_rrr::sum(unsigned long long i) {
  if (i == 0 ) return 0;
  unsigned long long index = rv_select_1(i);
  return rv_rank_0(index);
}

unsigned long long compressed_vector_rrr::search(unsigned long long j) {
  unsigned long long pos = rv_select_0(j);
  return rv_rank_1(pos);
}

size_t compressed_vector_rrr::get_size_in_bytes() {
  return sdsl::size_in_bytes(rv) + sdsl::size_in_bytes(rv_rank_0) +
                                   sdsl::size_in_bytes(rv_rank_1) +
                                   sdsl::size_in_bytes(rv_select_0) +
                                   sdsl::size_in_bytes(rv_select_1);
}

unsigned long long compressed_vector_rrr::get_count() {
  return rv_rank_1(rv.size()) + 1;
}

unsigned long long compressed_vector_rrr::get_next_64(unsigned long long j) {
  // std::string binary = std::bitset<64>(rv.get_int(j)).to_string(); //to binary
  // std::cout<<binary<<"\n";
  return rv.get_int(j);
}




compressed_vector_sd::compressed_vector_sd(vector<unsigned long long> &v) {

  // unsigned long long sm = 0;
  // for (auto it = v.begin(); it != v.end(); ++it)
  //   sm += (*it);
  // sdsl::bit_vector bv(sm + v.size() - 1, 0);
  //
  // sm = 0;
  // for (unsigned long long i=0; i<v.size()-1; ++i) {
  //   sm += v[i];
  //   bv[i+sm] = 1;
  // }

  rv = sdsl::sd_vector<>(v.begin(), v.end());

  rv_rank_0 = sdsl::sd_vector<>::rank_0_type(&rv);
  rv_rank_1 = sdsl::sd_vector<>::rank_1_type(&rv);
  rv_select_0 = sdsl::sd_vector<>::select_0_type(&rv);
  rv_select_1 = sdsl::sd_vector<>::select_1_type(&rv);

}


unsigned long long compressed_vector_sd::sum(unsigned long long i) {
  if (i == 0 ) return 0;
  unsigned long long index = rv_select_1(i);
  return rv_rank_0(index);
}

unsigned long long compressed_vector_sd::search(unsigned long long j) {
  unsigned long long pos = rv_select_0(j);
  return rv_rank_1(pos);
}

size_t compressed_vector_sd::get_size_in_bytes() {
  return sdsl::size_in_bytes(rv) + sdsl::size_in_bytes(rv_rank_0) +
                                   sdsl::size_in_bytes(rv_rank_1) +
                                   sdsl::size_in_bytes(rv_select_0) +
                                   sdsl::size_in_bytes(rv_select_1);
}

unsigned long long compressed_vector_sd::get_count() {
  return rv_rank_1(rv.size()) + 1;
}
