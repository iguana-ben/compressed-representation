#include "partition.h"

#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <set>
#include <functional>

#include "vec_tools.h"
#include "entropy_utils.h"

// TODO:
// -sprawdz indeksy w do_partition zero


// float h0_splitted(const vector<vector<int>> &sps) {
//
// 	std::map<vector<int>, int> mm;
// 	for(auto it = sps.begin(); it != sps.end(); ++it) {
// 		mm[*it]++;
// 	}
//
// 	float size = 0.0f;
// 	for(auto it = mm.begin(); it != mm.end(); ++it ) {
// 		size += float((it->second))*log2( float(sps.size())/float(it->second) );
// 	}
// 	return size;
//
// }
//
//
// //  i-factor start, l - length
// double quasi_optimal_partiton::factor_cost_zero(int i, int l, int max_order, vector< vector<int> > &freqs) {
//
// 	double cost = 0.0;
//
// 	for (int j=0; j < l; ++j) {
// 		int context_len = std::min(max_order, j);
//
// 		int start_index = i + j - context_len;
// 		if (context_len == 0) {
// 			cost += log2(double(freqs[context_len][i+j])/double(freqs[context_len+1][i+j]));
// 		} else {
// 			cost += log2(double(freqs[context_len][start_index])/double(freqs[context_len+1][start_index]));
// 		}
// 	}
//
// 	return cost;
//
// }
//
//
// //  i-factor start, lprev - previous factor length, l - length
// double quasi_optimal_partiton::factor_cost_first(
// 	int i, int lprev, int l, int max_order, vector< vector<int> > &freqs) {
//
// 	double cost = 0.0;
//
//
// 	if(i - lprev < 0 )
// 		throw new string("cannot compute cost: out of range\n");
//
// 	cost = factor_cost_zero(i-lprev, lprev + l, max_order, freqs);
// 	cost -= factor_cost_zero(i-lprev, lprev, max_order, freqs);
//
// 	return cost;
//
// }
//
// // 0 1 2 3 4 |
// //
// // i = 0 j = 4 context_len = 2
// vector< vector<int> > quasi_optimal_partiton::do_partition(
// 	const vector<int>& vec, const int max_factor_length,
// 	int enable_length_penalty, int dict_limit) {
//
// 	// const int max_factor_length = 5;
// 	const int order = max_factor_length-1;
//
// 	vector< vector<int> > result;
//
// 	set<vector<int>> top_set;
// 	auto freqs = fill_freq(vec, max_factor_length + 1, dict_limit, top_set);
// 	// for(auto it = top_set.begin(); it != top_set.end(); ++it) {
// 	// 	for(unsigned int j = 0; j < (*it).size(); ++j)
// 	// 		std::cout<<char((*it)[j]);
// 	// 	std::cout<<" | "<<(*it).size()<<" ";
// 	// 	std::cout<<"\n";
// 	// }
// 	// std::cout<<"dict end\n";
// 	// std::cout<<"ok so far\n";
//
// 	bool enable_set = (dict_limit != -1);
//
// 	double infty = 1000.0*vec.size(); //it is enough infty > vec.size()*log(alphabet size)
// 	vector<double> cost(vec.size() + 3, infty);
// 	vector<int> best_choice(vec.size() + 3, 0);
// 	cost[0] = 0; //0 - empty string
//
// 	for (int i=0; i < int(vec.size()) ; ++i) {
//
// 		cost[i+1] = infty*2.0;
// 		best_choice[i] = 1;
//
// 		for (int f = 1; f <=std::min(i+1, max_factor_length); ++f) {
//
// 			// if w is not in the most frequent set then so wa is not
// 			if (enable_set && f > 1) {
// 				if ( top_set.find(subvec(vec,i+1-f, i)) == top_set.end() )
// 					break;
// 			}
//
// 			double cc = cost[i+1-f] + factor_cost_zero(i+1-f, f, order, freqs);
//
// 			if (enable_length_penalty) {
// 				// double length_penalty = log2(double(max_factor_length)/f);
// 				// double length_penalty = f*log2(double(max_factor_length)/f);
// 				double length_penalty = 2.0*log2(double(f));
// 				cc += length_penalty;
// 			}
//
// 			if(cc < cost[i+1]) {
// 				cost[i+1] = cc;
// 				best_choice[i] = f;
// 			}
// 		}
//
// 	}
// 	int index = vec.size()-1;
// 	while(index > 0) {
// 		result.push_back(subvec(vec, index-best_choice[index]+1, index));
// 		index = index-best_choice[index];
// 	}
//
// 	std::reverse(result.begin(), result.end());
//
// 	return result;
// }
//
//
// vector< vector<int> > quasi_optimal_partiton::do_partition_first(
// 	const vector<int>& vec, const int max_factor_length,
// 	int enable_length_penalty, int dict_limit) {
//
// 	//TODO: co jesli wektor jest za krótki, a order za duży?
//
// 	// int dict_limit = 50000;
//
// 	const int order = std::min(int(max_factor_length*2 - 1), int(vec.size() - 1));
//
// 	vector< vector<int> > result;
//
// 	set<vector<int>> top_set;
// 	auto freqs = fill_freq(vec, order+1, dict_limit, top_set);
//
// 	double infty = 1000.0*vec.size(); //it is enough that infty > vec.size()*log(alphabet size)
// 	vector<vector<double>> cost(vec.size() + 3);
// 	vector<vector<int>> best_choice(vec.size() + 3);
//
// 	bool enable_set = (dict_limit != -1);
//
// 	for (unsigned int i = 0; i < cost.size(); ++i ) {
// 		cost[i] = vector<double>(order+1, infty);
// 		best_choice[i] = vector<int>(order+1, -1);
// 	}
//
// 	cost[0][1] = 0;
// 	best_choice[0][0] = 0;
//
// 	for (int i=1; i < int(vec.size()); ++i) {
//
// 		for (int factor_len = 1; factor_len <= max_factor_length; ++ factor_len) {
//
// 			if (enable_set && factor_len > 1) {
// 				if ( top_set.find(subvec(vec, i+1-factor_len, i)) == top_set.end() )
// 					break;
// 			}
//
// 			double best_cost = infty*vec.size();
// 			int b_choice = -1;
// 			// int ind = i-1;
// 			for (int prev_len = 1; prev_len <= max_factor_length; ++prev_len) {
// 				if (i - prev_len - factor_len + 1 >= 0 && i < int(vec.size()) ) {
// 					double cc = cost[i - factor_len][prev_len] +
// 						factor_cost_first(i - factor_len + 1, prev_len, factor_len, order, freqs);
//
// 					if (enable_length_penalty) {
// 						double length_penalty = log2(double(max_factor_length)/factor_len);
// 						cc += length_penalty;
// 					}
//
// 					if(cc < best_cost) {
// 						best_cost = cc;
// 						b_choice = prev_len;
// 					}
// 				}
// 			}
//
// 			cost[i][factor_len] = best_cost;
// 			best_choice[i][factor_len] = b_choice;
// 		}
// 	}
//
// 	double min_cost = infty*vec.size()*2.0;
// 	int c_len = -1;
// 	for (int i=1; i<= max_factor_length; ++i) {
// 		if(cost[vec.size()-1][i] <= min_cost) {
// 			min_cost = cost[vec.size()-1][i];
// 			c_len = i;
// 		}
// 	}
//
//
//
// 	double sanity_check = 0.0;
//
// 	int index = vec.size()-1;
// 	while(index > 0) {
// 		int prev_len = best_choice[index][c_len];
// 		sanity_check += factor_cost_first(index-c_len+1, prev_len, c_len, order, freqs);
//
// 		result.push_back(subvec(vec, index-c_len+1, index));
// 		int temp = c_len;
// 		c_len = best_choice[index][c_len];
// 		index -= temp;
// 	}
//
// 	result.push_back(subvec(vec, 0, 0));
//
// 	std::reverse(result.begin(), result.end());
//
// 	std::cout<<" Total cost: " << (long long) min_cost <<std::endl;
//
// 	return result;
// }
//
// // result[i][j] -- # of substrings equal to vec[i...i+j-1] in vec
// vector< vector<int> > quasi_optimal_partiton::fill_freq(
// 		const vector<int>& vec, int l, int topk, set<vector<int>> &top_set) {
//
// 	vector< vector<int> > result;
// 	vector<int> res; //res[i][0] = vec.size()
// 	for (unsigned int i = 0; i < vec.size(); ++i)
// 		res.push_back(vec.size());
//
// 	result.push_back(res);
//
// 	std::map<vector<int>, int> total_freqs;
//
// 	for (int i = 0; i < l; ++i) {
// 		vector<int> res;
// 		std::map< vector<int>, int> mm;
//
// 		// a bit hacky, works similarly to union-find
// 		// if x == y then freq[y] == freq[x]
// 		// i think order matters ( 0 to vec.size())
// 		auto cmp = [i, &vec] (const int& x, const int& y) {
// 			for (int z = 0; z <= i; ++z) {
// 				if (vec[x+z] != vec[y+z])
// 					return vec[x+z] < vec[y+z];
// 			}
// 			return false;
// 		};
// 		std::map<int, int, std::function<bool(const int&, const int&)> > freq_map(cmp);
//
// 		for (int j=0; j < int(vec.size())-i; ++j) {
// 			// mm[subvec(vec, j, j+i)]++;
// 			freq_map[j]++;
//
// 			if (i != 0 && i != l - 1 && topk != -1)
// 				total_freqs[subvec(vec, j, j+i)]++;
// 			// std::cout<<mm[subvec(vec, j, j+i)]<<"\n";
// 		}
//
// 		// for (int j=0;j<i; ++i)
// 			// result.push_back(-1);
//
// 		for (int j = 0; j < int(vec.size())-i; ++j) {
// 			// res.push_back(mm[subvec(vec, j, j+i)]);
// 			res.push_back(freq_map[j]);
//
// 			// std::cout<<mm[subvec(vec, j, j+i)] << " "<<freq_map[j]<<"\n";
// 		}
//
// 		result.push_back(res);
// 		std::cout<<i<<"\n";
// 	}
//
// 	vector<std::pair<int, vector<int>>> top_vector;
// 	for ( auto it = total_freqs.begin(); it != total_freqs.end(); ++it) {
// 		top_vector.push_back(make_pair(it->second, it->first));
// 	}
// 	std::sort(top_vector.begin(), top_vector.end(), std::greater<std::pair<int, vector<int>>>());
// 	for (int i=0; i < topk && i < int(top_vector.size()) ; ++i)
// 		top_set.insert(top_vector[i].second);
//
// 	return result;
// }


vector< string > quasi_optimal_partiton_string::best_shift_partition(
		const string& vec, int factor_length) {
	vector<string> result[factor_length];

	int best_index = 0;
	double best_cost;
	for (int i=0; i<factor_length; ++i) {

		int j = i;
		if (i != 0)
			result[i].push_back( vec.substr(0, i) );

		while ( j < int(vec.size()) ) {
			result[i].push_back( subvec(vec, j, std::min(j+factor_length-1,int(vec.size()-1) )) );
			j += factor_length;
		}

		if (i == 0) {
			best_cost = h0_splitted(result[i]);
			best_index = 0;
		} else {
			if (h0_splitted(result[i]) < best_cost) {
				best_index = i;
				best_cost = h0_splitted(result[i]);
			}
		}
	}
	return result[best_index];
}

// vector< vector<int> > quasi_optimal_partiton::naive_partition(
// 		const vector<int>& vec, int factor_length) {
// 	vector< vector<int> > result;
//
// 	int j = 0;
// 	while ( j < int(vec.size()) ) {
// 		result.push_back( subvec(vec, j, std::min(j+factor_length-1, int(vec.size()-1) )) );
// 		j += factor_length;
// 	}
//
// 	return result;
// }
//
//
// std::pair<std::map<int, vector<int>>, vector<int>> quasi_optimal_partiton::replace_partition(
// 		const vector<vector<int>> &part) {
// 	std::map<int, vector<int>> dictionary_map;
// 	std::map<vector<int>, int> r_dictionary_map;
// 	vector<int> text;
// 	int letter = 0;
// 	for(int i=0; i < int(part.size()); ++i) {
// 		if (r_dictionary_map.find(part[i]) == r_dictionary_map.end()) {
// 			text.push_back(letter);
// 			r_dictionary_map[part[i]] = letter;
// 			dictionary_map[letter] = part[i];
// 			letter++;
// 		} else {
// 			text.push_back(r_dictionary_map[part[i]]);
// 		}
// 	}
//
// 	return make_pair(dictionary_map, text);
// }
//
// double quasi_optimal_partiton::zero_order_partition_size(const vector<vector<int>> &part) {
// 	auto pp = replace_partition(part);
// 	//print_vec(pp.second); std::cout<< "\n";
// 	//std::cout<<part.size()<<" a\n";
// 	int sm =0;
// 	for (auto it = part.begin(); it != part.end(); ++it){
// 		// print_vec(it->second); std::cout<<"\n"
// 		sm += it->size();
// 		// std::cout<< it->second.size() << "a\n";
// 	}
// 	std::cout<<pp.first.size()<<" "<<pp.second.size()<<" "<<float(sm)/float(part.size())<<" distinct/count/avg_len\n";
//
// 	std::set<vector<int> > dict;
// 	for (auto it = part.begin(); it != part.end(); ++it )
// 		dict.insert(*it);
//
// 	for (auto it = dict.begin(); it != dict.end(); ++it ){
// 		// std::cout<<to_str_from_vec(*it)<<"\n";
// 	}
// 	// for ()
//
// 	return h0_splitted(part);
// }
//
// double quasi_optimal_partiton::first_order_partition_size(const vector<vector<int>> &part) {
// 	auto pp = replace_partition(part);
// 	//print_vec(pp.second); std::cout<< "\n";
// 	//std::cout<<part.size()<<" a\n";
// 	int sm =0;
// 	for (auto it = part.begin(); it != part.end(); ++it){
// 		// print_vec(it->second); std::cout<<"\n"
// 		sm += it->size();
// 		// std::cout<< it->second.size() << "a\n";
// 	}
// 	std::cout<<pp.first.size()<<" "<<pp.second.size()<<" "<<float(sm)/float(part.size())<<" distinct/count/avg_len\n";
//
// 	// std::set<vector<int> > dict;
// 	// for (auto it = part.begin(); it != part.end(); ++it )
// 	// 	dict.insert(*it);
//
// 	// for (auto it = dict.begin(); it != dict.end(); ++it ){
// 	// 	// std::cout<<to_str_from_vec(*it)<<"\n";
// 	// }
// 	// // for ()
// 	return bit_size(1, pp.second);
// 	// return h0_splitted(part);
// }
