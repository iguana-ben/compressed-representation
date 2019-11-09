#ifndef COMPRESSED_VECTOR_RRR_H
#define COMPRESSED_VECTOR_RRR_H

#include <vector>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>


using std::vector;

/*
codes vector v = {a_1, a_2, ..., a_n}
as 0^{a_1} 1 0^{a_2} 1 ... 1 0^{a_n}
*/
class compressed_vector_rrr {
  sdsl::rrr_vector<> rv;
  sdsl::rrr_vector<>::rank_0_type rv_rank_0;
  sdsl::rrr_vector<>::rank_1_type rv_rank_1;
  sdsl::rrr_vector<>::select_0_type rv_select_0;
  sdsl::rrr_vector<>::select_1_type rv_select_1;
public:
  compressed_vector_rrr() {} ;
  compressed_vector_rrr(const compressed_vector_rrr &cv) {
    rv = cv.rv;
    rv_rank_0 = sdsl::rrr_vector<>::rank_0_type(&rv);
    rv_rank_1 = sdsl::rrr_vector<>::rank_1_type(&rv);
    rv_select_0 = sdsl::rrr_vector<>::select_0_type(&rv);
    rv_select_1 = sdsl::rrr_vector<>::select_1_type(&rv);
  }

  compressed_vector_rrr(vector<unsigned long long> &v);
  unsigned long long sum(unsigned long long i);
  unsigned long long search(unsigned long long j);
  unsigned long long get_next_64(unsigned long long j); //v[j] must be 1
  size_t get_size_in_bytes();
  unsigned long long get_count();

};

class compressed_vector_sd {
  sdsl::sd_vector<> rv;
  sdsl::sd_vector<>::rank_0_type rv_rank_0;
  sdsl::sd_vector<>::rank_1_type rv_rank_1;
  sdsl::sd_vector<>::select_0_type rv_select_0;
  sdsl::sd_vector<>::select_1_type rv_select_1;
public:
  compressed_vector_sd() {} ;
  compressed_vector_sd(const compressed_vector_sd &cv) {
    rv = cv.rv;
    rv_rank_0 = cv.rv_rank_0;
    rv_rank_1 = cv.rv_rank_1;
    rv_select_0 = cv.rv_select_0;
    rv_select_1 = cv.rv_select_1;
  }

  compressed_vector_sd(vector<unsigned long long> &v);
  unsigned long long sum(unsigned long long i);
  unsigned long long search(unsigned long long j);
  size_t get_size_in_bytes();
  unsigned long long get_count();

};


#endif
