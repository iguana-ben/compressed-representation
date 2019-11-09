#include "entropy_utils.h"

#include <cmath>
#include <map>
#include <unordered_map>

#include <set>

#include <iostream>

#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

using std::map;
using std::pair;
using std::set;

double bit_size(int order, const vector<int> &str)
{
  vector<int> str2;

  for(int i=0; i<order; ++i)
  	str2.push_back(-1);

  str2.insert(str2.end(), str.begin(), str.end());

  // gp_hash_table<pair<int, vector<int> >, int> contx_ppb;
  // gp_hash_table<vector<int>, int> contx_sums;
  map<pair<int, vector<int> >, int > contx_ppb;
  map<vector<int>, int> contx_sums;

  vector<int> current_context;

  for(int i=0; i<order; ++i)
  	current_context.push_back(str2[i]);

  for(int i=order; i < int(str2.size()); ++i){

  	contx_sums[current_context]++;
  	contx_ppb[make_pair(str2[i], current_context) ]++;

  	current_context.push_back(str2[i]);
  	current_context.erase(current_context.begin());
  }

  double sum = 0.0;

  for(auto b = contx_ppb.begin(); b != contx_ppb.end(); b++) {
  	// sum += contx_sums
  	double pl = 0.0;
  	double lg = 0.0;
  	if(contx_sums[b->first.second] != b->second)
      lg = log2( (double(contx_sums[b->first.second]) + pl )/double(b->second) );

  	// print_vec(b->first.second);cout<<"\n";
  	// cout<<b->first.first<<" "<<b->second<<" "<<contx_sums[b->first.second]<<"\n";
  	// cout<<double(b->second)/double(contx_sums[b->first.second])<<"|\n";
  	// sum += contx_sums.begin()
  	sum += lg*double(b->second);
  }
  // cout<<"\ncontext_count: "<<contx_ppb.size()<<" "<< contx_ppb.size() * log2( 54 ) + sum <<" \n";
  return sum;
}

double bit_size(int order, const string &str)
{
  gp_hash_table<string, long int > contx_ppb;
  gp_hash_table<string, long int> contx_sums;

  string current_context = "";

  for(int i=0; i<order; ++i)
  	current_context += str[i];

  for(int i=order; i < int(str.size()); ++i){

  	contx_sums[current_context]++;
    contx_ppb[current_context + str[i]]++;

  	current_context += str[i];
  	current_context.erase(current_context.begin());
  }

  double sum = 0.0;

  for(auto b = contx_ppb.begin(); b != contx_ppb.end(); b++) {
  	double lg = 0.0;
    string ctx = b->first;
    ctx.pop_back();
    long int cnt = contx_sums[ctx];
  	if(cnt != b->second)
      lg = log2( (double(cnt) )/double(b->second) );
  	sum += lg*double(b->second);
  }
  return sum;
}

double h0_vec(const vector<int> &str) {
  gp_hash_table<int, long int > symbol_cnt;
  for (unsigned long int i=0; i <str.size(); ++i)
    symbol_cnt[str[i]]++;

  double sum = 0.0;
  for(auto it = symbol_cnt.begin(); it != symbol_cnt.end(); it++) {
    long int cnt = it->second;
    double lg = 0.0;
    if (cnt != (long int)str.size())
      lg = log2( (double(str.size()))/double(cnt) );
    sum += lg*double(cnt);
  }

  return sum;
}

double h1_vec(const vector<int> &str) {

  gp_hash_table<int, long int > symbol_cnt;
  std::map<pair<int, int>, long int > pair_cnt;

  for (unsigned long int i=1; i <str.size(); ++i) {
    symbol_cnt[str[i-1]]++;
    pair_cnt[std::make_pair(str[i-1], str[i])]++;
  }

  double sum = 0.0;
  for(auto it = pair_cnt.begin(); it != pair_cnt.end(); it++) {
    double lg = 0.0;
    int context = it->first.first;
    long int context_cnt = symbol_cnt[context];
    if (context_cnt != it->second)
      lg = log2( (double(context_cnt))/double(it->second) );
    sum += lg*double(it->second);
  }

  return sum;
}

int count_fo_dictionary_sizes(const vector<int> &str)
{
  map<int, set<int> > sizes;
  int sum = 0;

  for (unsigned int i=0; i<str.size()-1; ++i) {
    sizes[str[i]].insert(str[i+1]);
  }
  int count_single = 0;
  for (auto it = sizes.begin(); it != sizes.end(); ++it) {
    sum += it->second.size();
    count_single += (it->second.size() == 1);
    // std::cout << it->first<<" "<< it->second.size() << "\n";
  }
  std::cout<<"sum, single: "<<sum<<" "<<count_single<<"\n";
  // cout<<"\ncontext_count: "<<contx_ppb.size()<<" "<< contx_ppb.size() * log2( 54 ) + sum <<" \n";
  return sum;
}


unsigned long long count_fo_contexts(const vector<int> &str)
{
  set<pair<int, int>> contexts;
  for (unsigned long long i = 0; i < str.size()-1; ++i )
    contexts.insert(std::make_pair(str[i], str[i+1]));
  return contexts.size();
}
