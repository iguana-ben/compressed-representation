#ifndef VEC_TOOLS
#define VEC_TOOLS

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

using std::vector;
using std::string;


vector<int> subvec (const vector<int>& vec, int i, int j);
vector<int> to_vec(std::string str);

const int max_size = 1024*1024*1024; //1gb...
vector<int> load_from_file(const std::string filename);

// template<typename T> 
//void print_vec(const vector<int> &vv);

template<typename T> 
void print_vec(const vector<T> &vv) {
	for(int i=0; i<vv.size(); ++i)
		std::cout<<vv[i]<<" ";
	// cout<<"\n";
}

string to_str_from_vec(const vector<int> &vv);

#endif