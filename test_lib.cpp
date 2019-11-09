#include "test_lib.h"

#include "partition.h"
#include "vec_tools.h"
#include "entropy_utils.h"
#include "algo_utils.h"
#include "entropy_coding.h"
#include "compressed_array.h"
#include "compressor.h"

#define DEBUG_TEST 0

void test::test_directory(string input_directory, string output_directory) {
  return;
}

// compare_stat_record test::test_file(string input_filename, int order) {
//     compare_stat_record result;
//     result.entropies = std::vector<double>();
//     result.mean_entropies = std::vector<double>();
//
//     auto file = load_from_file(input_filename);
//     // for (int i=0; i <= order; ++i)
//
//     return result;
// }

partition_stat_record test::test_file(string input_filename, int max_factor_length,
    bool do_best_shift, bool do_zero, bool do_first, int enable_length_penalty,
    bool write_output, string output_filename, int dict_limit) {

  partition_stat_record record;

  auto file = load_from_file(input_filename);

  auto part_naive = quasi_optimal_partiton::naive_partition(file, max_factor_length);
  auto replaced_naive = quasi_optimal_partiton::replace_partition(part_naive);

  record.naive_zeroth_order_entropy =  bit_size(0, replaced_naive.second);
  record.naive_zeroth_order_seq_length = replaced_naive.second.size();
  record.naive_zeroth_order_dictionary_size = replaced_naive.first.size();

  record.naive_first_order_entropy =  bit_size(1, replaced_naive.second);
  record.naive_first_order_seq_length = replaced_naive.second.size();
  record.naive_first_order_dictionary_size = replaced_naive.first.size();

  // count_fo_dictionary_sizes(replaced_naive.second);

  if (do_best_shift) {
    auto part =
      quasi_optimal_partiton::best_shift_partition(file, max_factor_length);
    record.best_shift_seq_length_zero = part.size();
    auto replaced_partition = quasi_optimal_partiton::replace_partition(part);
    record.best_shift_entropy_zero = bit_size(0, replaced_partition.second);
    record.best_shift_dictionary_size_zero = replaced_partition.first.size();
  }

  if (do_zero) {
    auto part =
      quasi_optimal_partiton::do_partition(file, max_factor_length, enable_length_penalty, dict_limit);
    record.zeroth_order_seq_length = part.size();
    auto replaced_partition = quasi_optimal_partiton::replace_partition(part);
    record.zeroth_order_entropy = bit_size(0, replaced_partition.second);
    record.zeroth_order_dictionary_size = replaced_partition.first.size();
  }

  if (do_first) {
    auto part =
      quasi_optimal_partiton::do_partition_first(file, max_factor_length, enable_length_penalty, dict_limit);
    record.first_order_seq_length = part.size();
    auto replaced_partition = quasi_optimal_partiton::replace_partition(part);
    record.first_order_entropy = bit_size(1, replaced_partition.second);
    record.first_order_dictionary_size = replaced_partition.first.size();
  }

  return record;
}

partition_stat_record test::test_file_string(string input_filename, int max_factor_length,
    bool do_best_shift, bool do_zero, bool do_first, int enable_length_penalty,
    bool write_output, string output_filename, int dict_limit) {

  partition_stat_record record;

  memset(&record, 0, sizeof(partition_stat_record));

  string file = "";
  auto file_vec = load_from_file(input_filename);
  for (unsigned int i=0; i<file_vec.size(); ++i) {
    file.push_back(char(file_vec[i]));
  }
  file_vec = std::vector<int>();

  // {
  //   auto part_naive = quasi_optimal_partiton::naive_partition(file_vec, max_factor_length);
  //   auto replaced_naive = quasi_optimal_partiton::replace_partition(part_naive);
  //
  //   record.naive_zeroth_order_entropy =  bit_size(0, replaced_naive.second);
  //   record.naive_zeroth_order_seq_length = replaced_naive.second.size();
  //   record.naive_zeroth_order_dictionary_size = replaced_naive.first.size();
  //
  //   record.naive_first_order_entropy =  bit_size(1, replaced_naive.second);
  //   record.naive_first_order_seq_length = replaced_naive.second.size();
  //   record.naive_first_order_dictionary_size = replaced_naive.first.size();
  // }

  if (do_best_shift) {
    auto part_zero = quasi_optimal_partiton::best_shift_partition(0, file, max_factor_length);
    record.best_shift_seq_length_zero = part_zero.size();
    auto replaced_partition_zero = quasi_optimal_partiton::replace_partition(part_zero);
    record.best_shift_entropy_zero = h0_vec(replaced_partition_zero.second);
    record.best_shift_dictionary_size_zero = replaced_partition_zero.first.size();
    record.best_shift_dictionary_size_sum_zero = 0;
    for (auto it = replaced_partition_zero.first.begin(); it != replaced_partition_zero.first.end(); ++it)
      record.best_shift_dictionary_size_sum_zero += it->second.length();

    if (DEBUG_TEST) {
      std::cout<<"Checking best shift zero: ";
      bool is_equal = true;
      long int index = 0;
      for (uint i=0; i< part_zero.size(); ++i)
        for (uint j=0; j < part_zero[i].length(); ++j)
          is_equal = is_equal && part_zero[i][j] == file[index++];

      if (is_equal)
        std::cout<<"OK\n";
      else
        std::cout<<"Not equal\n";
    }

  }

  if (do_best_shift) {
    auto part_first = quasi_optimal_partiton::best_shift_partition(1, file, max_factor_length);
    record.best_shift_seq_length_first = part_first.size();
    auto replaced_partition_first = quasi_optimal_partiton::replace_partition(part_first);
    record.best_shift_entropy_first = h1_vec(replaced_partition_first.second);
    record.best_shift_dictionary_size_first = replaced_partition_first.first.size();
    record.best_shift_dictionary_size_sum_first = 0;
    for (auto it = replaced_partition_first.first.begin(); it != replaced_partition_first.first.end(); ++it)
      record.best_shift_dictionary_size_sum_first += it->second.length();

    record.best_shift_contexts = count_fo_contexts(replaced_partition_first.second);

    if (DEBUG_TEST) {
      std::cout<<"Checking best shift first: ";
      bool is_equal = true;
      long int index = 0;
      for (uint i=0; i< part_first.size(); ++i)
        for (uint j=0; j < part_first[i].length(); ++j)
          is_equal = is_equal && part_first[i][j] == file[index++];

      if (is_equal)
        std::cout<<"OK\n";
      else
        std::cout<<"Not equal\n";
    }

  }

  if (do_zero) {
    auto part =
      quasi_optimal_partiton::do_partition(file, max_factor_length, enable_length_penalty, dict_limit);
    record.zeroth_order_seq_length = part.size();
    auto replaced_partition = quasi_optimal_partiton::replace_partition(part);
    record.zeroth_order_entropy = h0_vec(replaced_partition.second);
    record.zeroth_order_dictionary_size = replaced_partition.first.size();
    record.zeroth_order_dictionary_size_sum = 0;
    for (auto it = replaced_partition.first.begin(); it != replaced_partition.first.end(); ++it)
      record.zeroth_order_dictionary_size_sum += it->second.length();

    if (DEBUG_TEST) {
      std::cout<<"Checking part zero: ";
      bool is_equal = true;
      long int index = 0;
      for (uint i=0; i< part.size(); ++i)
        for (uint j=0; j < part[i].length(); ++j)
          is_equal = is_equal && part[i][j] == file[index++];

      if (is_equal)
        std::cout<<"OK\n";
      else
        std::cout<<"Not equal\n";
    }

  }

  if (do_first ) {
    auto part =
      quasi_optimal_partiton::do_partition_first(file, max_factor_length, enable_length_penalty, dict_limit);
    record.first_order_seq_length = part.size();
    auto replaced_partition = quasi_optimal_partiton::replace_partition(part);
    record.first_order_entropy = h1_vec(replaced_partition.second);
    record.first_order_dictionary_size = replaced_partition.first.size();
    record.first_order_dictionary_size_sum = 0;
    for (auto it = replaced_partition.first.begin(); it != replaced_partition.first.end(); ++it)
      record.first_order_dictionary_size_sum += it->second.length();

    record.first_order_contexts = count_fo_contexts(replaced_partition.second);

    if (DEBUG_TEST) {
      std::cout<<"Checking part first: ";
      bool is_equal = true;
      long int index = 0;
      for (uint i=0; i< part.size(); ++i)
        for (uint j=0; j < part[i].length(); ++j)
          is_equal = is_equal && part[i][j] == file[index++];

      if (is_equal)
        std::cout<<"OK\n";
      else
        std::cout<<"Not equal\n";
    }

  }

  return record;

}

void test::test_structure_zeroth(string input_filename, int max_factor_length,
    int enable_length_penalty, int dict_limit) {
    auto file = load_from_file(input_filename);
    string file_str = "";
    for (uint i=0; i<file.size(); ++i)
      file_str.push_back(char(file[i]));


    char* buffer_uncompressed = new char[file.size() + 5];
    char* buffer_naive = new char[file.size() + 5];
    char* buffer_alg = new char[file.size() + 5];

    double bps_result[3];
    double random_result[3];
    double block_result[3];


    unsigned int total_to_read = 1000000;
    unsigned int loops = 1;
    vector<unsigned int> indices;

    srand((unsigned)time(NULL));

    for (unsigned int i=0; i<total_to_read; ++i)
      indices.push_back(rand() % (file.size()-5) );

    unsigned int total_to_read_blocks = 1000;
    unsigned int block_length = 51200;
    unsigned int loops_block = 1;
    vector<unsigned int> indices_block;

    for (unsigned int i=0; i<total_to_read_blocks; ++i)
      indices_block.push_back(rand() % (file.size()-5-block_length) );


    auto compressed_file_naive = compressed_array_zeroth(file_str, max_factor_length, enable_length_penalty, 0);
    auto compressed_file_alg = compressed_array_zeroth(file_str, max_factor_length, enable_length_penalty, 1);

    bps_result[0] = 8.0;
    bps_result[1] = double(compressed_file_naive.get_size_in_bits())/file.size();
    bps_result[2] = double(compressed_file_alg.get_size_in_bits())/file.size();

    std::cout<<"starting time tests: \n";
    auto chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops; ++l) {
      for (unsigned int i=0; i<total_to_read; ++i) {
        buffer_uncompressed[i] = file_str[indices[i]];
         // ((char*)(file_str.data()))[indices[i]];
      }
    }
    auto chrono_end = std::chrono::high_resolution_clock::now();
    random_result[0] =  std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;

    chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops; ++l) {
      for (unsigned int i=0; i<total_to_read; ++i) {
        buffer_naive[i] = compressed_file_naive.get(indices[i]);;
      }
    }
    chrono_end = std::chrono::high_resolution_clock::now();
    random_result[1] =  std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;


    chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops; ++l) {
      for (unsigned int i=0; i<total_to_read; ++i) {
        buffer_alg[i] = compressed_file_alg.get(indices[i]);
      }
    }
    chrono_end = std::chrono::high_resolution_clock::now();
    random_result[2] =  std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;



    std::cout<<"bps: "<<bps_result[0] <<" " << bps_result[1] <<" " << bps_result[2] << std::endl;
    std::cout.precision(10);
    std::cout<<std::fixed<<"random time: "<<random_result[0] <<" " << random_result[1] <<" " << random_result[2] << std::endl;

    chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops_block; ++l) {
      for (unsigned int i=0; i<total_to_read_blocks; ++i) {
        memcpy(buffer_uncompressed, file_str.data() +  indices_block[i], block_length);
      }
    }
    chrono_end = std::chrono::high_resolution_clock::now();
    block_result[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;

    chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops_block; ++l) {
      for (unsigned int i=0; i<total_to_read_blocks; ++i) {
        compressed_file_naive.read_block(buffer_naive, indices_block[i], block_length);
      }
    }
    chrono_end = std::chrono::high_resolution_clock::now();
    block_result[1] = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;


    chrono_begin = std::chrono::high_resolution_clock::now();
    for (unsigned int l=0; l<loops_block; ++l) {
      for (unsigned int i=0; i<total_to_read_blocks; ++i) {
        compressed_file_alg.read_block(buffer_alg, indices_block[i], block_length);
      }
    }
    chrono_end = std::chrono::high_resolution_clock::now();
    block_result[2] = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end-chrono_begin).count() / 1000000000.0;

    std::cout<<std::fixed<<"block time: "<<block_result[0]<<" " <<block_result[1]<<" " <<block_result[2] << std::endl;

    // std::cout<<" check correctness: \n";
    // for (uint j=0; j<50; ++j)
    //   std::cout<<buffer_uncompressed[j] <<" " << buffer_naive[j] <<" " << buffer_alg[j] << "\n";
    //

    std::cout<<" check correctness: \n";
    for (uint j=0; j<50; ++j)
      std::cout<<buffer_uncompressed[j] <<" " << buffer_naive[j] <<" " << buffer_alg[j] << "\n";


   for (uint j=0; j<1000; ++j)
    std::cout<<file_str[j] <<" " << compressed_file_naive.get(j)
             << " " <<compressed_file_alg.get(j) <<"\n";

}

compare_stat_record test::entropy_compare(string input_filename, int m, bool do_fo) {
  compare_stat_record result;
  result.entropies = vector<double>();
  result.mean_entropies = vector<double>();
  result.mean_entropies_d  = vector<double>();

  auto file = load_from_file(input_filename);
  string file_str = "";
  for (uint i=0; i<file.size(); ++i)
    file_str.push_back(char(file[i]));

  double sm = 0.0;
  vector<double> prefix_sums_entropy;

  for (int i=0; i < m + m*do_fo; ++i) {
    double bs = bit_size(i, file_str);
    result.entropies.push_back(bs);
    sm += bs;
    result.mean_entropies.push_back(sm/(i+1));
    prefix_sums_entropy.push_back(sm);

    std::cout<<"Finished: "<<i<<std::endl;
  }

  if (do_fo) {
    for (int i=1; i <= m; ++i) {
      result.mean_entropies_d.push_back( (prefix_sums_entropy[2*i-1]-prefix_sums_entropy[i-1])/i );
    }
  }

  return result;

}


void test::compress_test(string input_filename, int max_factor_length) {
  auto file = load_from_file(input_filename);
  string file_str = "";
  for (uint i=0; i<file.size(); ++i)
    file_str.push_back(char(file[i]));

  for (int method = 0; method < 2; ++method) {

    compressor::bps_data bps;

    vector<bool> compressed_data =
      compressor::compress_zeroth(max_factor_length, method, file_str, &bps);

    std::cout<<"File: " + input_filename + " Bps (total/huffman/dict): "
             << bps.bps << " " << bps.bitstring_bps <<" "<< bps.dict_bps <<"\n";

    std::cout<<"finished compression"<<std::endl;

    string decompressed = compressor::decompress_zeroth(compressed_data);

    if (DEBUG_TEST) {
      if (file_str == decompressed)
        std::cout<<"Files equal, OK\n";
      else
        std::cout<<"Files NOT equal, FAILED\n";
    }
  }
}

void test::compress_test_first(string input_filename, int max_factor_length) {
  auto file = load_from_file(input_filename);
  string file_str = "";
  for (uint i=0; i<file.size(); ++i)
    file_str.push_back(char(file[i]));

  for (int method = 0; method < 2; ++method) {
    compressor::bps_data bps;
    vector<bool> compressed_data =
      compressor::compress_first(max_factor_length, method, file_str, &bps);

    std::cout<<"File: " + input_filename + " Bps (total/huffman/dict): "
             << bps.bps << " " << bps.bitstring_bps <<" "<< bps.dict_bps <<"\n";


    std::cout<<"finished compression"<<std::endl;

    string decompressed = compressor::decompress_first(compressed_data);

    if (file_str == decompressed)
      std::cout<<"Files equal, OK\n";
    else
      std::cout<<"Files NOT equal, FAILED\n";
  }
}
