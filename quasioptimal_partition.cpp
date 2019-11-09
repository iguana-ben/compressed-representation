#include <iostream>
#include <string>

#include "partition.h"
#include "vec_tools.h"
#include "entropy_utils.h"
#include "test_lib.h"
#include "algo_utils.h"
#include "compressed_vector_rrr.h"
#include "entropy_coding.h"


using namespace std;



void test_structre(string name);
void test_compression(string name);
void test_stats(string name, bool first_or_second_test);

int main() {
	string name = "bible.txt";

	test_structre(name);
	test_compression(name);
	test_stats(name, 0);
	test_stats(name, 1);


	return 0;
}

void test_structre(string name) {
	test::test_structure_zeroth("test_files/" + name, 8, 1, -1);
}

void test_compression(string name) {
	for (int i=1; i<9; ++i) {
		std::cout<<"Testing H0, m: "<<i<<"\n";
		test::compress_test("test_files/" + name, i);
	}

	for (int i=1; i<5; ++i) {
		std::cout<<"Testing H1, m: "<<i<<"\n";
		test::compress_test_first("test_files/" + name, i);
	}

}

// used to generate entropy data (Table 1/2 in paper)
void test_stats(string name, bool first_or_second_test) {

	std::cout<<"File: "<<name<<std::endl;
	int m = 4;
	auto entropy = test::entropy_compare("test_files/" + name, m, first_or_second_test);
	for (uint i = 0; i < entropy.entropies.size(); ++i) {
			// std::cout<<"H"<<i<<" H_mean "<<(long int)(entropy.entropies[i])<<" "
							// <<(long int)(entropy.mean_entropies[i])<<" "<<(long int)(entropy.mean_entropies_d[i])<<std::endl;
			std::cout<<"H"<<i<<" H_mean "<<(long int)(entropy.entropies[i])<<" "
							<<(long int)(entropy.mean_entropies[i])<<" "<<std::endl;
	}
	for (uint i = 0; i<entropy.mean_entropies_d.size(); ++i) {
		std::cout<<"Mean2: "<<i+1<<" "<<(long int)(entropy.mean_entropies_d[i]) << std::endl;
	}
	// return 0;

	for (int x = 1; x <= m; ++x) {

	auto test_stat = test::test_file_string("test_files/" + name, x, 1, !first_or_second_test, first_or_second_test, 1, 0, "", -1);
	std::cout<<"Current max len:" << x << "\n";
	std::cout<<(long int)(test_stat.best_shift_entropy_zero)<<" "<<test_stat.best_shift_seq_length_zero<<" "
					 <<test_stat.best_shift_dictionary_size_zero <<" "<<test_stat.best_shift_dictionary_size_sum_zero<<" \n";
	std::cout<<(long int)(test_stat.zeroth_order_entropy)<<" "<<test_stat.zeroth_order_seq_length<<" "
					 <<test_stat.zeroth_order_dictionary_size<<" "<<test_stat.zeroth_order_dictionary_size_sum<<" \n";
	std::cout<<(long int)(test_stat.best_shift_entropy_first)<<" "<<test_stat.best_shift_seq_length_first
					 <<" "<<test_stat.best_shift_dictionary_size_first<<" "<<test_stat.best_shift_dictionary_size_sum_first
					 <<" "<<test_stat.best_shift_contexts<<" \n";
	std::cout<<(long int)(test_stat.first_order_entropy)<<" "<<test_stat.first_order_seq_length
					<<" "<<test_stat.first_order_dictionary_size<<" "<<test_stat.first_order_dictionary_size_sum
					<<" "<<test_stat.first_order_contexts<<" \n";
	}

}
