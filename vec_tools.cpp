#include "vec_tools.h"

#include <iomanip>

// vec[i..j]
vector<int> subvec (const vector<int>& vec, int i, int j) {

	return vector<int>(vec.begin() + i, vec.begin() + j + 1);
}

vector<int> to_vec(std::string str) {
	std::vector<int> res;
	for(int i=0; i < int(str.length()); ++i)
		res.push_back(str[i]);
	return res;
}


vector<int> load_from_file(const std::string filename) {

	std::ifstream fl(filename);

    fl.seekg( 0, std::ios::end );

    size_t len = fl.tellg();
    char *ret = new char[max_size];
    fl.seekg(0, std::ios::beg);

    fl.read(ret, std::min(max_size, int(len) ));

    vector<int> res;
    for(int i=0; i<std::min(max_size, int(len) ); ++i) {
    	res.push_back(ret[i]);
    }

    delete[] ret;

    return res;
}

string to_str_from_vec(const vector<int> &vv) {

		string str = "";
		for (int i=0; i < int(vv.size()) ; ++i) {
			str += char(vv[i]);
		}

		return str;
}
