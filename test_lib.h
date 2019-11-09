#ifndef TEST_LIB_H
#define TEST_LIB_H

#include <string>
#include <vector>

using std::string;

struct partition_stat_record {
	double naive_zeroth_order_entropy;
	int naive_zeroth_order_seq_length;
	int naive_zeroth_order_dictionary_size;

	double naive_first_order_entropy;
	int naive_first_order_seq_length;
	int naive_first_order_dictionary_size;

	double best_shift_entropy_zero;
	int best_shift_seq_length_zero;
	int best_shift_dictionary_size_zero;
	long int best_shift_dictionary_size_sum_zero;

	double best_shift_entropy_first;
	int best_shift_seq_length_first;
	int best_shift_dictionary_size_first;
	long int best_shift_dictionary_size_sum_first;
	long int best_shift_contexts;

	double zeroth_order_entropy;
	int zeroth_order_seq_length;
	int zeroth_order_dictionary_size;
	long int zeroth_order_dictionary_size_sum;

	double first_order_entropy;
	int first_order_seq_length;
	int first_order_dictionary_size;
	long int first_order_dictionary_size_sum;
	long int first_order_contexts;

};

struct compare_stat_record {
	std::vector<double> entropies; //first entropies
	std::vector<double> mean_entropies;
	std::vector<double> mean_entropies_d;

	int num_lz78_factors;
	double lz78_entropy;

	long int file_length;
};

class test {
public:
	static void test_directory(string input_directory, string output_directory);
	static partition_stat_record test_file(string input_filename,
			int max_factor_length, bool do_best_shift, bool do_zero, bool do_first,
			int enable_length_penalty, bool write_output, string output_filename, int dict_limit);
	static partition_stat_record test_file_string(string input_filename,
			int max_factor_length, bool do_best_shift, bool do_zero, bool do_first,
			int enable_length_penalty, bool write_output, string output_filename, int dict_limit);
	static void test_structure_zeroth(string input_filename, int max_factor_length,
	    int enable_length_penalty, int dict_limit);
	static compare_stat_record entropy_compare(string input_filename, int m, bool do_fo);
	static void compress_test(string input_filename, int max_factor_length);
	static void compress_test_first(string input_filename, int max_factor_length);
};


#endif
