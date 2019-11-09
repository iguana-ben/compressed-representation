#include "entropy_coding.h"
#include <cassert>

struct huffman_node
{
	huffman_node *left, *right;
	unsigned int letter;
	huffman_node(){left = NULL; right = NULL; }
};

void huffman_traverse(huffman_node* node, map<unsigned int, string> &res, string &path)
{
	if(node->left == NULL)
		res[node->letter] = path;
	else
	{
		path.push_back('0');
		huffman_traverse(node->left, res, path);
		path.pop_back();

		path.push_back('1');
		huffman_traverse(node->right, res, path);
		path.pop_back();
	}
}

void del_huffman_tree(huffman_node* node)
{
	if(node->left != NULL)
	{
		del_huffman_tree(node->left);
		del_huffman_tree(node->right);
	}
	delete node;
}


map<unsigned int, string> gen_huffman_codes( const vector<unsigned int> &text ) {

	map<unsigned int, string> dict;
	map<unsigned int, unsigned int> frequencies;
	for(unsigned int i=0; i < text.size(); ++i) {
		frequencies[text[i]]++;
	}

	std::priority_queue< pair<pair<unsigned int, int>, huffman_node*>,
						vector<pair<pair<unsigned int, int>, huffman_node*>>,
					std::greater< pair<pair<unsigned int, int>, huffman_node*>>> pq;

	int id = 0;
	map<unsigned int, unsigned int>::iterator it;
	for(it = frequencies.begin(); it != frequencies.end(); ++it)
	{
		// cout<<it->first<<" : "<<it->second<<"\n";
		huffman_node *node = new huffman_node();
		node->letter = it->first;
		pq.push( std::make_pair(std::make_pair(it->second, id++), node) );
	}

	while(pq.size() > 1)
	{
		auto n1 = pq.top(); pq.pop();
		auto n2 = pq.top(); pq.pop();
		huffman_node *new_node = new huffman_node();

		new_node->left = n1.second;
		new_node->right = n2.second;

		pq.push(std::make_pair(std::make_pair(n1.first.first + n2.first.first, id++), new_node));
	}

	if(frequencies.size() == 1)
		dict[text[0]] = "";
	else
	{
		string path = "";
		huffman_traverse(pq.top().second, dict, path);
	}
	return dict;
}

unsigned long long entropy_coding::code_to_binary(const string &code) {
	unsigned long long res = 0;
	unsigned long long mul = 1;
	string new_code = "1" + code;
	assert(new_code.length() <= 64); //todo add support for longer codes

	for (int i=int(new_code.length())-1; i>=0; --i) {
		if (new_code[i] == '1')
			res += mul;
		mul <<= 1;
	}

	return res;
}

string entropy_coding::binary_to_code(unsigned long long num) {
	string res = "";

	while(num) {
		res = char('0' + (num % 2)) + res;
		num >>= 1;
	}
	res.erase(res.begin());
	return res;
}

std::tuple<vector<unsigned long long>, vector<unsigned int>, map<unsigned int, string> >
	entropy_coding::get_huffman_codes_succinct(const vector<unsigned int> &text) {

  auto codes_map = gen_huffman_codes(text);
	vector<pair<unsigned long long, unsigned int> > codes_sorted;
	for (auto it = codes_map.begin(); it != codes_map.end(); ++it) {
		codes_sorted.push_back( std::make_pair( code_to_binary(it->second), it->first ));
		// std::cout<<it->second<<" "<<code_to_binary(it->second)<<" "<<binary_to_code(code_to_binary(it->second))<<"\n";
	}

	std::sort(codes_sorted.begin(), codes_sorted.end());

	auto result = make_tuple(vector<unsigned long long>(), vector<unsigned int>(), codes_map);
	for (unsigned int i = 0; i<codes_sorted.size(); ++i) {
		std::get<0>(result).push_back(codes_sorted[i].first);
		std::get<1>(result).push_back(codes_sorted[i].second);
		// std::cout <<codes_sorted[i].first<<" "<<codes_sorted[i].second<<"\n";
	}

	return result;

}

vector<std::tuple<vector<unsigned long long>, vector<unsigned int>, map<unsigned int, string> > >
	entropy_coding::get_huffman_codes_succinct_first_order(const vector<unsigned int> &text) {

	vector<std::tuple<vector<unsigned long long>,
				 vector<unsigned int>,
				 map<unsigned int, string> > > result;

	std::set<unsigned int> alphabet;
	for (uint i=0; i<text.size(); ++i)
		alphabet.insert(text[i]);

	unsigned int alphabet_size = alphabet.size();

	vector<vector<unsigned int> > Tsigmas(alphabet_size);

	for (uint i=1; i<text.size(); ++i) {
		Tsigmas[text[i-1]].push_back(text[i]);
	}

	for (uint i=0; i<alphabet_size; ++i) {
		//if last character is unique it would have empty string, we create placeholder, consumes few more bits, easier implementation
		if(Tsigmas[i].size() == 0)
			Tsigmas[i].push_back(0);

		result.push_back(get_huffman_codes_succinct(Tsigmas[i]));
	}

	return result;

}
