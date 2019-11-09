#ifndef ELIAS_CODING
#define ELIAS_CODING

class elias_coding {
	public:
		static vector<bool> encode(unsigned long long x) {

			vector<bool> bin;
			while(x > 0) {
				if(x % 2)
					bin.push_back(1);
				else
					bin.push_back(0);

				x/=2;
			}

			vector<bool> res = bin;

			res[bin.size()-1] = 1;
			for(uint i=1; i<bin.size(); ++i) {
				res[i-1] = 0;
				res.push_back(bin[bin.size()-1-i]);
			}
			return res;
		}

		static void decode(vector<bool> &vec, int &pos, unsigned long long &decoded_x) {
			int len = 0;
			while(!vec[pos+len])
				len++;
			// cout<<len<<"\n";
			unsigned long long decoded = 0;
			unsigned long long p2 = 1;
			// 0 1 2 3 4 5 6
			// 0 0 0 1 0 0
			// pos = 1 + 4
			for (int i=pos + 2*len; i >= pos + len; --i) {
				if(vec[i])
					decoded += p2;
				p2 *= 2;
			}

			decoded_x = decoded;
			pos += 2*len + 1;
		}
};

#endif
