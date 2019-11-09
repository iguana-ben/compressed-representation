#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <ext/pb_ds/assoc_container.hpp>


using namespace __gnu_pbds;

using std::vector;
using std::set;
using std::string;


// vector<int> subvec (const vector<int>& vec, int i, int j);

class quasi_optimal_partiton {
private:
	static double factor_cost_zero(int i, int l, int max_order, vector< vector<int> > &freqs);
	static double factor_cost_zero_string(long long i, long long l, long long max_order, const string &s, gp_hash_table<string, long long> &freqs);
	static double factor_cost_first(int i, int lprev, int l, int max_order, vector< vector<int> > &freqs);
	static double factor_cost_first_string(long long i, long long lprev, long long l, long long max_order, const string &s, gp_hash_table<string, long long> &freqs);
	static vector< vector<int> > fill_freq(const vector<int>& vec, int l, int topk, set<vector<int>> &top_set);
	// static gp_hash_table<string, int> fill_freq_string(const string& vec, int l, int topk, set<string> &top_set);
	static gp_hash_table<string, long long> fill_freq_string( const string& vec, long long l, long long topk, set<string> &top_set);
	static double calc_length_penalty(int f, int max_factor_length);

public:
	static std::pair<std::map<int, vector<int>>, vector<int>> replace_partition(const vector<vector<int>> &part);
	static std::pair<std::map<int, string>, vector<int>> replace_partition(const vector<string> &part);
	static double zero_order_partition_size(const vector<vector<int>> &part);
	static double first_order_partition_size(const vector<vector<int>> &part);
	static vector< vector<int> > naive_partition(const vector<int>& vec, int factor_length);
	static vector< vector<int> > best_shift_partition(const vector<int>& vec, int factor_length);
	static vector< vector<int> > do_partition(const vector<int>& vec, const int max_factor_length,
			int enable_length_penalty, int dict_limit);
	static vector< vector<int> > do_partition_first(const vector<int>& vec, const int max_factor_length,
			int enable_length_penalty, int dict_limit);

	//functions on strings are much faster but support alphabets to up to 256
	static vector< string > best_shift_partition(bool zero_or_one, const string& vec, int factor_length);
	static vector< string > do_partition(const string& vec, const int max_factor_length,
			int enable_length_penalty, int dict_limit);
	static vector< string > do_partition_first(const string& vec, const int max_factor_length,
			int enable_length_penalty, int dict_limit);

};


#endif
