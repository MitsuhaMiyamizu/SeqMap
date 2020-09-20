/*
probe.h - Header file for probe manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2007 -
*/

#ifndef PROBE_H
///Define this macro to prevent from including this header file more than once.
#define PROBE_H

#include "stl.h"
#include "type.h"
#include "string_operation.h"

const char DNA_C[4] = {'A', 'C', 'G', 'T'};
const int MAX_PROBE_LEN = 32;
int probe_len = 25;

inline string tostring(const int64 key, int len = probe_len){
	int i;
	int64 temp = key;
	char buf[MAX_PROBE_LEN + 1];
	for (i = 0; i < len; i++) {
		int last_char = (int)(temp & 3);
		buf[len - i - 1] = DNA_C[last_char];
		temp >>= 2;
	}
	buf[len] = 0;
	return string(buf);
}

inline int64 tonumber(const string &probe){
	int i;
	int64 result = 0;
	int probe_len = (int)probe.length();
	for (i = 0; i < probe_len; i++) {
		result <<= 2;
		char ch = probe[i];
		int j;
		for (j = 0; j < 4; j++) {
			if (ch == DNA_C[j]) break;
		}
		if (j == 4) {
			return -1;
		}
		result |= j;
	}
	return result;
}

inline char random_NT() {
	return (DNA_C[rand() % 4]);
}

inline bool fix_N_probe(string &probe, int max_N = 0) {
	int i;
	int probe_len = (int)probe.length();
	int num_fixed = 0;
	bool has_N = false;
	for (i = 0; i < probe_len; i++) {
		char ch = probe[i];
		int j;
		for (j = 0; j < 4; j++) {
			if (ch == DNA_C[j]) break;
		}
		if (j == 4) {
			has_N = true;
			if (num_fixed < max_N) {
				num_fixed++;
    			if (ch == 'N' || ch == 'n' || ch == '.') probe[i] = random_NT();
			}
		}
	}
	return has_N;
}

inline void fix_short_probe(string &probe, int max_N = 0) {
	int len = (int)probe.length();
	if (len >= probe_len) return;
	while ((int)probe.length() < probe_len) probe = probe + "N";
	if (probe_len - len <= max_N) fix_N_probe(probe, max_N);
}

inline void local_alignment(string &probe1, string &probe2, int num_insdel) {
	probe1 = toupper(probe1);
	probe2 = toupper(probe2);

	int i, j;
	int len1 = (int)probe1.length(), len2 = (int)probe2.length();
	string seq1, seq2;

	if (num_insdel == 0) {
		__ASSERT(len1 == len2, "internal error, len1 != len2.\n");
		for (i = 0; i < len1; i++) {
			if (probe1[i] == probe2[i]) {
				seq1.push_back(probe1[i]);
				seq2.push_back(probe2[i]);
			} else {
				seq1.push_back((char)tolower(probe1[i]));
				seq2.push_back((char)tolower(probe2[i]));
			}
		}
	} else {
		vector<vector<int> > d;
		int dim = 2 * max(len1, len2);
		d.resize(dim);
		for (i = 0; i < (int)d.size(); i++) d[i].resize(dim);
		for (i = 0; i <= len1; i++) d[i][0] = i;
		for (j = 0; j <= len2; j++) d[0][j] = j;

		for (i = 1; i <= len1; i++) {
			for (j = 1; j <= len2; j++) {
				int cost = (probe1[i-1] == probe2[j-1]) ? 0 : 1;
				d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1);
				d[i][j] = min(d[i][j], d[i-1][j-1]+cost);
			}
		}
		i = len1;
		j = len2;
		while (i > 0 || j > 0) {
			if (i > 0 && j > 0 && d[i][j] == d[i-1][j-1] + 1) {
				i--;
				j--;
				seq1.insert(seq1.begin(), (char)tolower(probe1[i]));
				seq2.insert(seq2.begin(), (char)tolower(probe2[j]));
			} else if (i > 0 && d[i][j] == d[i-1][j] + 1) {
				i--; 
				seq1.insert(seq1.begin(), (char)tolower(probe1[i]));
				seq2.insert(seq2.begin(), '_');
			} else if (j > 0 && d[i][j] == d[i][j-1] + 1) {
				j--;
				seq1.insert(seq1.begin(), '_');
				seq2.insert(seq2.begin(), (char)tolower(probe2[j]));
			} else if (i > 0 && j > 0 && d[i][j] == d[i-1][j-1]) {
				i--;
				j--;
				seq1.insert(seq1.begin(), probe1[i]);
				seq2.insert(seq2.begin(), probe2[j]);
			} else {
				panic("internal error: _edit_dist_pos\n");
			}
		}
	}
	probe1 = seq1;
	probe2 = seq2;
}

inline int _edit_dist3(const string &probe1, const string &probe2, int num_mismatch, int num_insdel, int &len) {
	int i, j;
	vector<vector<int> > d;
	int dim = max((int)probe1.length(), (int)probe2.length()) + num_insdel + 1;
	d.resize(dim);
	for (i = 0; i < (int)d.size(); i++) d[i].resize(dim);

	for (i = 0; i <= num_insdel; i++) d[i][0] = d[0][i] = i;
	int probe1_len = (int)probe1.length(), probe2_len = (int)probe2.length();
	for (i = num_insdel + 1; i <= probe1_len; i++) d[i][i - num_insdel - 1] = num_mismatch + 1;
	for (i = num_insdel + 1; i <= probe2_len; i++) d[i - num_insdel - 1][i] = num_mismatch + 1;

	int min_d = probe2_len;
	int min_d_i = 0;
	for (i = 1; i <= probe1_len; i++) {
		j = (i <= num_insdel) ? 1 : i - num_insdel;
		while (j <= i + num_insdel && j <= probe2_len) {
			int cost = (probe1[i-1] == probe2[j-1]) ? 0 : 1;
			d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1);
			d[i][j] = min(d[i][j], d[i-1][j-1]+cost);
			if (j == probe2_len && d[i][j] < min_d) {
				min_d = d[i][j];
				min_d_i = i;
			}
			if (i == j && d[i][j] > num_mismatch + num_insdel) return (num_mismatch + 1);
			j++;
		}
	}
	if (min_d == d[probe2_len][probe2_len]) len = probe2_len;
	else len = min_d_i;
	return min_d;
}

inline int _edit_dist2(const string &probe1, const string &probe2, int num_mismatch, int num_insdel, int &len) {
	int d[MAX_PROBE_LEN*2][MAX_PROBE_LEN+1];
	int i, j;
	for (i = 0; i <= num_insdel; i++) d[i][0] = d[0][i] = i;
	for (i = num_insdel + 1; i <= len; i++) d[i][i - num_insdel - 1] = num_mismatch + 1;
	for (i = num_insdel + 1; i <= probe_len; i++) d[i - num_insdel - 1][i] = num_mismatch + 1;

	int min_d = probe_len;
	int min_d_i = 0;
	for (i = 1; i <= len; i++) {
		j = (i <= num_insdel) ? 1 : i - num_insdel;
		while (j <= i + num_insdel && j <= probe_len) {
			int cost = (probe1[i-1] == probe2[j-1]) ? 0 : 1;
			d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1);
			d[i][j] = min(d[i][j], d[i-1][j-1]+cost);
			if (j == probe_len && d[i][j] < min_d) {
				min_d = d[i][j];
				min_d_i = i;
			}
			if (i == j && d[i][j] > num_mismatch + num_insdel) return (num_mismatch + 1);
			j++;
		}
	}
	if (min_d == d[probe_len][probe_len]) len = probe_len; 
	else len = min_d_i;
	return min_d;
}

inline int _edit_dist(const string &probe1, const string &probe2) {
	int d[MAX_PROBE_LEN+1][MAX_PROBE_LEN+1];
	int i, j;
	for (i = 0; i <= probe_len; i++) d[i][0] = d[0][i] = i;

	for (i = 1; i <= probe_len; i++) {
		for (j = 1; j <= probe_len; j++) {
			int cost = (probe1[i-1] == probe2[j-1]) ? 0 : 1;
			d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1);
			d[i][j] = min(d[i][j], d[i-1][j-1]+cost);
		}
	}
	return d[probe_len][probe_len];
}

inline vector<uchar> _edit_dist_pos(const string &probe1, const string &probe2) {
	int d[MAX_PROBE_LEN*2][MAX_PROBE_LEN*2];
	int i, j;
	int len1 = (int)probe1.length(), len2 = (int)probe2.length();
	for (i = 0; i <= len1; i++) d[i][0] = i;
	for (j = 0; j <= len2; j++) d[0][j] = j;

	for (i = 1; i <= len1; i++) {
		for (j = 1; j <= len2; j++) {
			int cost = (probe1[i-1] == probe2[j-1]) ? 0 : 1;
			d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1);
			d[i][j] = min(d[i][j], d[i-1][j-1]+cost);
		}
	}
	vector<uchar> result;
	result.clear();
	i = len1;
	j = len2;
	while (i > 0 || j > 0) {
		if (i > 0 && j > 0 && d[i][j] == d[i-1][j-1] + 1) {
			i--;
			j--;
			result.insert(result.begin(), (uchar)(j));
		} else if (i > 0 && d[i][j] == d[i-1][j] + 1) {
			i--;
			result.insert(result.begin(), (uchar)(j));
		} else if (j > 0 && d[i][j] == d[i][j-1] + 1) {
			j--;
			result.insert(result.begin(), (uchar)(j));
		} else if (i > 0 && j > 0 && d[i][j] == d[i-1][j-1]) {
			i--;
			j--;
		} else {
			panic("internal error: _edit_dist_pos\n");
		}
	}
	return result;
}

inline int _edit_dist2(const int64 key1, const int64 key2, int num_mismatch, int num_insdel, int &len) {
	return _edit_dist2(tostring(key1, len), tostring(key2), num_mismatch, num_insdel, len);
}

inline int _edit_dist(const int64 key1, const int64 key2) {
	return _edit_dist(tostring(key1), tostring(key2));
}

inline vector<uchar> _edit_dist_pos(const int64 key1, const int64 key2, int len=probe_len) {
	return _edit_dist_pos(tostring(key1, len), tostring(key2));
}

inline int64 get_reverse(const int64 key, int len = probe_len){
	int64 result = 0;
	int64 temp = (~key)&((((int64)(1))<<(len*2))-1);
	int i;
	for (i = 0; i < len; i++) {
		result <<= 2;
		result |= (temp & 3);
		temp >>= 2;
	}
	return result;
}

inline string get_reverse(const string &probe){
	int probe_len = (int)probe.length();
	return (tostring(get_reverse(tonumber(probe), probe_len), probe_len));
}

inline string get_reverse2(const string &probe){
	string result = "";
	string temp = toupper(probe);
	for (int i = (int)probe.length() - 1; i >= 0; i--) {
		if (probe[i] == 'A') result += "T";
		else if (probe[i] == 'T') result += "A";
		else if (probe[i] == 'C') result += "G";
		else if (probe[i] == 'G') result += "C";
		else result += probe[i];
	}
	return result;
}

static int mismatch_table[65536];

inline void build_mismatch_table() {
	int i;
	for (i = 0; i < 65536; i++) {
		int result = 0;
		int temp = i;
		while (temp) {
			if (temp & 3) result++;
			temp >>= 2;
		}
		mismatch_table[i] = result;
	}
}

inline int _mismatch2(const string &seq1, const string &seq2, int num_mismatch) {
	int result = 0; 
	for (int i = 0; i < min((int)seq1.length(), (int)seq2.length()); i++) {
		if (seq1[i] != seq2[i]) {
			result++;
			if (result > num_mismatch) return result;
		}
	}
	if ((int)seq1.length() > (int)seq2.length()) result += (int)seq1.length() - (int)seq2.length();
	else result += (int)seq2.length() - (int)seq1.length();
	return result;
}

inline int _mismatch2(const int64 key1, const int64 key2, int num_mismatch) {
	int64 temp = key1 ^ key2;
	int result = mismatch_table[temp & 65535];
	if (result > num_mismatch) return result;
	result += mismatch_table[(temp >> 16) & 65535];
	if (result > num_mismatch) return result;
	result += mismatch_table[(temp >> 32) & 65535];
	if (result > num_mismatch) return result;
	result += mismatch_table[(temp >> 48) & 65535];
	return result;
}

inline int _mismatch(const int64 key1, const int64 key2) {
	int64 temp = key1 ^ key2;
	int result = 0;
	int i;
	for (i = 0; i < probe_len; i++) {
		if ((temp & 3) != 0) result++;
		temp >>= 2;
	}
	return result;
}

inline vector<uchar> _mismatch_pos(const int64 key1, const int64 key2) {
	int64 temp = key1 ^ key2;
	int i;
	vector<uchar> result;
	result.clear();
	for (i = 0; i < probe_len; i++) {
		if ((temp & 3) != 0) result.insert(result.begin(), (uchar)(probe_len - 1 - (i)));
		temp >>= 2;
	}
	return result;
}

inline void correct_key(int64 &key, int len = probe_len) {
	key &= (((int64)(1) << (len * 2)) - 1);
}

inline int mismatch(const int64 key1, const int64 key2, int num_mismatch, int num_insdel, int &len){
	if (num_insdel > 0) return _edit_dist2(key1, key2, num_mismatch, num_insdel, len);
	else return _mismatch2(key1, key2, num_mismatch);
}

inline int mismatch(const string &seq1, const string &seq2, int num_mismatch, int num_insdel, int &len) {
	if (num_insdel > 0) return _edit_dist3(seq1, seq2, num_mismatch, num_insdel, len);
	else return _mismatch2(seq1, seq2, num_mismatch);
}

inline vector<uchar> mismatch_pos(const int64 key1, const int64 key2, int num_insdel, int &len){
	if (num_insdel > 0) return _edit_dist_pos(key1, key2, len);
	else {
		__ASSERT(len == probe_len, "internal error: wrong probe len.\n");
		return _mismatch_pos(key1, key2);
	}
}

#endif //PROBE_H
