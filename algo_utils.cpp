#include "algo_utils.h"

#include <set>
#include <iostream>

std::vector<std::vector<int>> lz78_factor(std::vector<int> &v) {

	std::set< std::vector<int>> factors;

	std::vector<std::vector<int> > result;
	std::vector<int> first;
	first.push_back(v[0]);
	result.push_back(first);
	for (unsigned int i=1; i < v.size();) {
		std::vector<int> factor;
		while (i < v.size() ) {
			factor.push_back(v[i++]);
			if (factors.find(factor) == factors.end() ) {
				factors.insert(factor);
				result.push_back(factor);
				break;
			}
		}
	}

	// std::cout<<result.size()<< " \n";
	// for (unsigned int i=0; i < result.size(); ++i ) {
	// 	for (unsigned int j=0; j<result[i].size(); ++j) {
	// 		std::cout<<result[i][j];
	// 	}
	// 	std::cout<<"|";
	// }
	// std::cout<<"\n";


	return result;
}
