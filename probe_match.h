/*
probe_match.h - Header file for probe matching
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
May, 2007 -
*/

#ifndef PROBE_MATCH_H
///Define this macro to prevent from including this header file more than once.
#define PROBE_MATCH_H

#include "stl.h"
#include "probe.h"
#include "system_utils.h"
#include "timer.h"
#include "string_operation.h"
#include "file_operation.h"
#include "fasta.h"

const int max_mismatch = 5;
const int max_insdel = 5;

#ifdef __DEBUG
	static int64 num_total_search_step = 0;
	static int64 num_search = 0;
	static int64 num_max_search_step = 0;

	static int64 total_jump1 = 0, total_jump2 = 0;
	static int64 max_jump1 = 0, max_jump2 = 0;
	static int64 total_insertion = 0, total_detection = 0;
#endif

class search_result_struct{
public:
	int trans_id;
	int trans_coord;
	int64 trans_key;
	int trans_key_len;
	string trans_seq;
	int probe_index;
	int num_mismatch;
	bool reverse_strand;

	bool reported;
};

typedef void search_result_handler(vector<search_result_struct> &results);

class probe_list;

class probe_info{
public:
	string id;
	string seq;
	int64 key;
	bool valid;
	bool rematch;
	int repeat; //-1 unique, >=0 identical to another probe
};

bool operator<(const probe_info& a, const probe_info& b) {
    return a.key < b.key;
}

class probe_match{
public:
	bool verbose;
	bool exact_mismatch;
	int num_mismatch;
	bool search_repeats;
	bool search_both_strands;
	bool duplicate_probes;
	bool cluster_identical_probes;

	int cut_start, cut_end;
	bool match_shorter_probes;
	bool match_N;

	bool filter_results;
	bool filter_selected_probes;
	string filter_selected_probes_filename;
	bool output_filtered_probes;
	bool filter_low_quality_probes;
//for hash
	bool use_hash;
	bool hash_1bp_mismatch;
//for probe lists
	int limit_num_part;
//for probe_list
	bool fast_index;
	double fast_index_fraction;
	bool interpolation_search;
	bool store_key;
	bool shift_mask;
//for insdel
	int num_insdel;
//for output
	bool zero_indexed;

	string probe_input_file_name;

	size_t used_memory;
	size_t available_memory;
	size_t extra_memory_per_probe;
	int num_probes;
	int num_valid_probes;
	int nn;
	uint hash_magic;
	uint* hash_table;
	vector<probe_info> probes;
	vector<string> transcripts;
	vector<probe_list*> probe_lists;
	probe_match();
	~probe_match();
	bool initialize();
	bool initialize(const string &fa_file_name);
	bool initialize(vector<string> &probe_sequences, vector<string> &probe_tags);
	bool initialize(vector<string> &probe_sequences) {
		vector<string> probe_tags;
		int n = (int)probe_sequences.size();
		probe_tags.resize(n);
		for (int i = 0; i < n; i++) probe_tags[i] = int2str(i + 1);
		return initialize(probe_sequences, probe_tags);		
	}
	bool search_key(const string &key, vector<uint> &indexes) { return search_key(tonumber(key), indexes);}
	bool search_key(const int64 key, vector<uint> &indexes);
	void search_results_prehandler(vector<search_result_struct> &results, search_result_handler handler);
	bool search_file(const string &fa_file_name, search_result_handler handler);
//for hash
	inline int64 get_probe(const uint probe_index);
	inline bool hash_insert_probe(const uint probe_index, const uint pos, const uint changeto);
	inline uint hash_detect_probe(const int64 key, vector<uint> &indexes);
};

class probe_list{
public:
	probe_match* pm;
	vector<int64> masks;
	int64 mask0_table[65536][4];
	int mask0_shift0, mask0_shift1, mask0_shift2;
	int fast_index_shift;
	int size_fast_index;
	int* fast_indexes;
	vector<int64> key_list;
	int* index_list;

	probe_list(){	
		fast_index_shift = 0;
		masks.clear();
		index_list = NULL;
		fast_indexes = NULL;
	}

	~probe_list(){
		if (NULL != index_list) delete[] index_list;
		if (NULL != fast_indexes) delete[] fast_indexes;
	}

	inline int64 shift_key(const int64 key, const int64 mask) {
		int64 result = 0;
		int64 m = mask;
		int64 k = key;
		int64 pos = 1;
		while (m) {
			result |= (k & m & 1) * pos;
			pos <<= (m & 1);
			m >>= 1;
			k >>= 1;
		}
		return result;
	}

	inline int64 shift_mask0_key(const int64 key) {
		int64 result = 0;
		result += mask0_table[key & 65535][0];
		result += (mask0_table[(key >> 16) & 65535][1] << mask0_shift0);
		result += (mask0_table[(key >> 32) & 65535][2] << mask0_shift1);
		result += (mask0_table[(key >> 48) & 65535][3] << mask0_shift2);
		return result;
	}

	inline void build_mask0_tables(){
		int i, j;
		for (i = 0; i < 65536; i++) {
			for (j = 0; j < 4; j++) {
				int64 submask = ((masks[0] >> (j * 16)) & 65535);
				mask0_table[i][j] = shift_key(i, submask);
			}
		}
		mask0_shift0 = count1(masks[0] & 65535);
		mask0_shift1 = mask0_shift0 + count1((masks[0] >> 16) & 65535);
		mask0_shift2 = mask0_shift1 + count1((masks[0] >> 32) & 65535);
	}

	inline int count1(int64 key) {
		int64 k = key;
		int result = 0;

		while (k) {
			result += int(k&1);
			k >>= 1;
		}
		return result;
	}

	inline int64 mask0_key(const int64 key) {
		if (pm->shift_mask) return shift_mask0_key(key);
		else return masks[0] & key;
	}

	inline int64 mask_key(const int64 key, const int64 mask) {
		if (pm->shift_mask) return (shift_key(key, mask)); 
		else return mask & key;
	}

	inline void set_mask(vector<int64> &_masks){
		masks = _masks;
		build_mask0_tables();
		if ((int64)(masks.size()) * pm->num_probes > ((int64)(1) << 31)) panic("internal error: out of integer bound.");
		int i, j;
		vector<pair<int64, int> > temp_probe_list;
		for (i = 0; i < pm->num_probes; i++)
			for (j = 0; j < (int)masks.size(); j++) {
				if (pm->probes[i].valid) temp_probe_list.push_back(pair<int64, int>(mask_key(pm->probes[i].key, masks[j]), j * pm->num_probes + i));
			}
		__ASSERT((int)temp_probe_list.size() == num_items(), "internal error: wrong num_items.\n");
		sort(temp_probe_list.begin(), temp_probe_list.end());
		index_list = new int[num_items()];
		for (i = 0; i < num_items(); i++) {
			index_list[i] = temp_probe_list[i].second;
			if (pm->store_key) key_list.push_back(temp_probe_list[i].first);
		}
		temp_probe_list.clear();
		if (pm->fast_index) {
			int64 temp = shift_key(masks[0], masks[0]);
			temp++;
			fast_index_shift = 0;
			while(temp > num_items() / pm->fast_index_fraction){
				temp >>= 1;
				fast_index_shift++;
			}
			size_fast_index = (int)(shift_key(masks[0], masks[0]) >> fast_index_shift) + 1;
			fast_indexes = new int[size_fast_index];
			fast_indexes[0] = 0;
			int current = 1;
			for (i = 0; i < num_items(); i++) {
				int temp;
				if (masks.size() > 1) {
					__ASSERT(pm->num_insdel > 0 && pm->shift_mask, "internal error: masks.size > 0.\n");
				}
				if (pm->shift_mask) temp = (int)(get_key(i) >> fast_index_shift);
				else temp = (int)(shift_key(pm->probes[get_probe_index(i)].key, masks[0]) >> fast_index_shift);
				if (temp >= current) {
					for (j = current; j <= temp; j++) {
						fast_indexes[j] = i;
					}
					current = temp+1;
				}
			}
			for (j = current; j < size_fast_index; j++) fast_indexes[j] = num_items() - 1;
		}
	}

	inline int64 get_key(const int index) {
		if (pm->store_key) return key_list[index];
		else return (masks.size() > 1) ? mask_key(get_probe(index), get_mask(index)) : mask_key(pm->probes[index_list[index]].key, masks[0]);
	}

	inline int64 get_probe(const int index) {
		return pm->probes[get_probe_index(index)].key;
	}

	inline int64 get_mask(const int index) {
		return masks[get_mask_index(index)];
	}

	inline int get_probe_index(const int index) {
		return (masks.size() > 1) ? index_list[index] % pm->num_probes : index_list[index];
	}

	inline int get_mask_index(const int index) {
		if (masks.size() > 1) {
			return index_list[index] / pm->num_probes;
		} else return 0;
	}

	inline int num_items() {
		return pm->num_valid_probes * (int)masks.size();
	}

	inline bool search_probe(const int64 probe, vector<uint> &indexes){
		int left, right, mid;
		int64 subprobe = (probe >> (pm->num_insdel * 2));
		int64 key = mask0_key(subprobe);
		
		if (pm->fast_index) {
			int temp = int(shift_mask0_key(subprobe) >> fast_index_shift);
			left = fast_indexes[temp];
			if (temp < size_fast_index - 1) right = fast_indexes[temp + 1];
			else right = num_items() - 1;
		} else {
			left = 0;
			right = num_items() - 1;
		}
		
		mid = left;
#ifdef __DEBUG
		int64 count_search_step = 0;
#endif
		while((get_key(left)) < key && (get_key(right)) >= key) {
			__ASSERT((left < right), "internal error in search");
#ifdef __DEBUG
			count_search_step++;
#endif
			if (pm->interpolation_search) {
				mid = left + (int)((((double)key - get_key(left)) * (right - left)) / (get_key(right) - get_key(left)));
				__ASSERT((left <= mid && mid <= right), "internal error in search");
			} else mid = (left + right) / 2;
			if ((get_key(mid)) > key) 
				right = mid - 1;
			else if ((get_key(mid)) < key)
				left = mid + 1;
			else {
				left = mid;
				break;
			}
		}
#ifdef __DEBUG
		num_total_search_step += count_search_step;
		num_search++;
		if (count_search_step > num_max_search_step) num_max_search_step = count_search_step;
#endif
		if ((get_key(left)) != key) 
			return false;
		else {
			while (left > 0 && (get_key(left-1)) == key)
				left--;
			while (left < num_items() && (get_key(left)) == key) {
				int len = probe_len + pm->num_insdel;
				int score = mismatch(probe, get_probe(left), pm->num_mismatch, pm->num_insdel, len);

				if (pm->exact_mismatch && score == pm->num_mismatch) indexes.push_back(get_probe_index(left));
				else if (score <= pm->num_mismatch) indexes.push_back(get_probe_index(left));

				left++;
			}
			return true;
		}
	}
};

inline uint joaat_hash(const uchar *key, const uint len)
 {
	 uint hash = 0;
     uint i;

     for (i = 0; i < len; i++) {
         hash += key[i];
         hash += (hash << 10);
         hash ^= (hash >> 6);
     }
     hash += (hash << 3);
     hash ^= (hash >> 11);
     hash += (hash << 15);
     return hash;
 }

inline int64 change_probe(const int64 probe, const uint pos, const uint changeto){
	return((probe & (((int64)(1) << (probe_len * 2)) - 1 - ((int64)(3) << (pos * 2)))) | ((int64)(changeto) << (pos * 2)));
}

inline uint encode_probe_index(const uint probe_index, const uint pos, const uint changeto){
	return (probe_index | (pos << 25) | (changeto << 30));
}

inline uint decode_probe_index(const uint probe_index){
	return (probe_index & ((1 << 25) - 1));
}

inline int64 probe_match::get_probe(const uint probe_index){
	uint index = decode_probe_index(probe_index);
	uint pos = (probe_index >> 25) & 31;
	uint changeto = (probe_index >> 30) & 3;
	__ASSERT(index < (uint)num_probes,"get probe error");
	__ASSERT(pos < (uint)probe_len, "get probe error");
	return (change_probe(probes[index].key, pos, changeto));
}

inline bool probe_match::hash_insert_probe(const uint probe_index, const uint pos, const uint changeto) {
	uint encoded_probe_index = encode_probe_index(probe_index, pos, changeto);
	int64 temp1 = get_probe(encoded_probe_index);
	uint index = joaat_hash((uchar*)(&temp1), 8) % nn;
	uint index_inc = 1;
	uint origin = index;
#ifdef __DEBUG
	int64 jump = 1;
#endif
	while (hash_table[index] != (uint)(-1)) {
#ifdef __DEBUG
		jump++;
#endif
		index += index_inc;
		index %= nn;
		if (origin == index)
			return false;
	}
#ifdef __DEBUG
	if (jump > max_jump1) max_jump1 = jump;
	total_jump1 += jump;
	total_insertion++;
#endif
	hash_table[index] = encoded_probe_index;
	return true;
}

inline uint probe_match::hash_detect_probe(const int64 key, vector<uint> &indexes) {
	uint index = joaat_hash((uchar*)(&key), 8) % nn;
	uint index_inc = 1;
	uint origin = index;
	indexes.clear();
#ifdef __DEBUG
	int64 jump = 1;
#endif
	while (hash_table[index] != (uint)(-1)) {
#ifdef __DEBUG
		jump++;
#endif
		int64 temp;
		temp = get_probe(hash_table[index]);

		if (temp == key) {
			indexes.push_back(hash_table[index]);
		}
		index += index_inc;
		index %= nn;
		if (origin == index) panic("circle in hash table");
	}
#ifdef __DEBUG
	if (jump > max_jump2) max_jump2 = jump;
	total_jump2 += jump;
	total_detection++;
#endif
	if (indexes.size() == 0) 
		return (uint)(-1);
	else
		return indexes[0];
}

inline probe_match::probe_match(){
	verbose = true;
	limit_num_part = 0;
	exact_mismatch = false;
	search_repeats = true;
	search_both_strands = true;
	duplicate_probes = true;
	cluster_identical_probes = true;

	cut_start = cut_end = 0;
	match_shorter_probes = false;
	match_N = true;

	filter_results = true;
	filter_selected_probes = false;
	filter_selected_probes_filename = "";
	output_filtered_probes = false;
	filter_low_quality_probes = false;

	use_hash = false;
	hash_1bp_mismatch = true;
	hash_table = NULL;
	num_mismatch = 0;
	num_probes = 0;
	num_valid_probes = 0;

	fast_index = true;
	fast_index_fraction = 4;
	interpolation_search = false;
	shift_mask = false;
	store_key = true;
	
	used_memory = 0;
	available_memory = 0;
	extra_memory_per_probe = 0;
	num_insdel = 0;

	zero_indexed = false;
}

inline probe_match::~probe_match(){
	int i;
	for (i = 0; i < (int)probe_lists.size(); i++)
		delete probe_lists[i];
	if (use_hash) {
		if (hash_table != NULL) {
			delete[] hash_table;
			hash_table = NULL;
		}
	}
}

inline static void _insdel_masks(vector<int64> &masks, int64 mask, int num_insdel){
	if (num_insdel == 0) return;
	int i;
	for (i = 0; i < probe_len; i++) {
		if ((mask & ((int64)(3) << (i * 2))) == 0) {
			int64 left = ((mask >> (i * 2)) << (i * 2));
			int64 right = mask - left;
			int64 new_mask = left + (right >> 2);
			if (find(masks.begin(), masks.end(), new_mask) == masks.end()) masks.push_back(new_mask);
			_insdel_masks(masks, new_mask, num_insdel-1);
			new_mask = left + (right << 2);
			if (find(masks.begin(), masks.end(), new_mask) == masks.end()) masks.push_back(new_mask);
			_insdel_masks(masks, new_mask, num_insdel-1);
		}
	}
}

inline bool probe_match::initialize(){
	int i, j, k;

	for (i = 0; i < (int)probe_lists.size(); i++)
		delete probe_lists[i];
	if (use_hash) {
		if (hash_table != NULL) {
			delete[] hash_table;
			hash_table = NULL;
		}
	}

	num_probes = (int)probes.size();
	if(num_probes <= 0) panic("no probe.\n");
	num_valid_probes = 0;
	for (i = 0; i < num_probes; i++) {
		if (probes[i].valid) num_valid_probes++;
	}
	if(num_valid_probes <= 0) panic("no valid probe.\n");
	if (verbose) cout << num_probes << " probes imported. " << num_valid_probes << " are valid." << endl;
	if (verbose) cout << "checking resources...";

	if (is32system()) {
		if (verbose) cout << "32-bit version." << endl;
	} else {
		if (verbose) cout << "64-bit version." << endl;
	}
	if (available_memory == 0) available_memory = get_available_memory_size();
	if (verbose) cout << "available memory: " << (int)(available_memory >> 20) << "MB.\n";

	build_mismatch_table();

	if (!search_both_strands) duplicate_probes = false;
	if (num_insdel > 0) {
		use_hash = false;
		shift_mask = true;
		if (search_both_strands) duplicate_probes = true;
	} else {
		filter_results = false;
	}

	if (num_mismatch > max_mismatch) panic("num_mismatch > max_mismatch.\n");
	if (num_insdel > max_insdel) panic("num_insdel > max_insdel.\n");
	if (num_insdel > num_mismatch) panic("num_insdel > num_mismatch.\n");

	if (duplicate_probes) {
		if (verbose) cout << "generating reversed probes for searching reverse strand\n";
		for (i = 0; i < num_probes; i++) {
			probe_info temp_probe_info = probes[i];
			if (temp_probe_info.valid) {
				if (temp_probe_info.rematch && (int)temp_probe_info.seq.length() > probe_len) {
					string probe = get_reverse2(temp_probe_info.seq);
					probe = probe.substr(0, probe_len);
					if (match_N) {
						fix_N_probe(probe, num_mismatch);
					}
					temp_probe_info.key = tonumber(probe);
					if (temp_probe_info.key == -1) {
						temp_probe_info.valid = false;
						probes[i].valid = false;
						num_valid_probes--;
					}
				} else {
					temp_probe_info.key = get_reverse(temp_probe_info.key);
				}
			}
			probes.push_back(temp_probe_info);
		}
		__ASSERT(num_probes * 2 == (int)probes.size(), "internal error: wrong num_probes.\n");
		num_probes = (int)probes.size();
		int old_num_valid_probes = num_valid_probes;
		num_valid_probes = 0;
		for (i = 0; i < num_probes; i++) {
			if (probes[i].valid) num_valid_probes++;
		}
		__ASSERT(old_num_valid_probes * 2 == num_valid_probes, "internal error: wrong num_valid_probes.\n");
		if (verbose) cout << "now " << num_probes << " probes. " << num_valid_probes << " are valid." << endl;
	}

	used_memory = probes.size() * (sizeof(probe_info) + probes[num_probes-1].id.length() + probes[num_probes-1].seq.length()); //probes
	used_memory += probes.size() * extra_memory_per_probe;

	if (use_hash) {
		if (verbose) cout << "building hash table...";
		double NN;
		if (num_mismatch > 0 && hash_1bp_mismatch)
			NN = num_probes * probe_len * 3 * 2.0;
		else
			NN = num_probes * 2.0;

		int logNN = (int)ceil(log(NN) / log(double(2.0)));
		if (logNN >= 32) panic("too many probes");
		nn = 1 << logNN;
		used_memory += (size_t)nn * sizeof(uint);
		if (verbose) cout << "estimated memory usage: " << (int)(used_memory >> 20) << "MB.\n";
		if (used_memory > available_memory) {
			panic("not enough memory.\n");
		}
		if (hash_table != NULL) {
			delete hash_table;
			hash_table = NULL;
		}
		hash_table = new uint[nn];
	rebuild:
		hash_magic = rand();
		for (i = 0; i < nn; i++) {
			hash_table[i] = uint(-1);
		}
		for (i = 0; i < num_probes; i++) {
			if (verbose && i % (num_probes / 20) == 0)  cout << '.'<< flush;
			int inserted_own = false;
			if (!probes[i].valid) continue;
			int64 key = probes[i].key;
			if (num_mismatch == 0 || !hash_1bp_mismatch) {
				if (!hash_insert_probe(i, 0, (uint)(key & 3))) {
					if (verbose) cout << " circle found, rebuild" << endl;
					goto rebuild;
				}
				continue;
			}
			for (j = 0; j < probe_len; j++) {
				for (k = 0; k < 4; k++) {
					int64 temp = get_probe(encode_probe_index(i, j, k));
					if (temp == key) {
						if (inserted_own) continue;
						inserted_own = true;
					}
					if (!hash_insert_probe(i, j, k)) {
						if (verbose) cout << " circle found, rebuild" << endl;
						goto rebuild;
					}
				}
			}
		}
		if (verbose) {
			if (verbose) cout << num_probes << " probes inserted, size:" << sizeof(uint) * nn << " bytes" << endl;
#ifdef __DEBUG
			if (verbose) cout << "average insertion time:" << (double)total_jump1 / total_insertion << endl;
			if (verbose) cout << "maximum insertion time:" << max_jump1 << endl;
#endif
			if (verbose) cout << "testing hash table...";
		}
		for (i = 0; i < num_probes; i++) {
			if (!probes[i].valid) continue;
			int64 temp = probes[i].key;
			vector<uint> indexes;
			bool in = false;
			if (hash_detect_probe(temp, indexes) != (uint)(-1)) {
				for (j = 0; j < (int)indexes.size(); j++)
					if (decode_probe_index(indexes[j]) == uint(i)) {
						in = true;
						break;
					}
			}
			if (!in) panic("error");
			if (num_mismatch > 0 && hash_1bp_mismatch) {
				uint j = rand()%probe_len;
				uint k = rand()%4;
				temp = change_probe(temp, j, k);
				vector<uint> indexes;
				bool in = false;
				if (hash_detect_probe(temp, indexes) != (uint)(-1)) {
					for (j = 0; j < (int)indexes.size(); j++)
						if (decode_probe_index(indexes[j]) == uint(i)) {
							in = true;
							break;
						}
				}
				if (!in) panic("error");
			}
			if (verbose && (i % (num_probes/20) == 0)) cout << "."<< flush;
		}
		if (verbose) {
			if (verbose) cout << "passed." << endl;
#ifdef __DEBUG
			if (verbose) cout << "average query time:" << (double)total_jump2 / total_detection << endl;
			if (verbose) cout << "maximum query time:" << max_jump2 << endl;
#endif
		}
	} else {
	//compute num_part
		int num_part = num_mismatch + 1;
		int effect_len = probe_len - num_insdel;
		while (true) {
			if (num_part == effect_len) break;
			int num_char_per_part = effect_len / num_part;
			int effect_char = num_char_per_part * (num_part - num_mismatch);
			if (effect_char*2 >= (int)(log((double)(num_probes))/log((double)2.0)))
				break;
			num_part++;
		}

		if(limit_num_part > 0 && num_part > limit_num_part) num_part = limit_num_part;

		if (verbose) cout << "Split probe into " << num_part << " parts." << endl;

		int num_lists = 1;
		for (i = 0; i < num_mismatch; i++)
			num_lists *= (num_part - i);
		for (i = 0; i < num_mismatch; i++)
			num_lists /= (i + 1);

		vector<vector<int64> > mask_sets;
		mask_sets.clear();

	//build probe lists
		int part_change[max_mismatch];

		for (j = 0; j < num_mismatch; j++)
			part_change[j] = j;
		
		while (true) {
			int64 m = 0;
			for (j = 0; j < num_part; j++) {
				for (k = 0; k < num_mismatch; k++)
					if (part_change[k] == j) break;
				if (k < num_mismatch) continue;
				int i1 = effect_len * j / num_part + num_insdel;
				int i2 = effect_len * (j + 1) / num_part + num_insdel;
				for (k = i1; k < i2; k++)
					m |= (((int64)(1) << (k*2)) | ((int64)(1) << (k*2+1)));
			}

			vector<int64> masks;
			masks.push_back(m);
			mask_sets.push_back(masks);

			int cur_pos = num_mismatch - 1;
			if (cur_pos == -1) break;
			part_change[cur_pos] = part_change[cur_pos] + 1;
			while (part_change[cur_pos] == num_part - num_mismatch + cur_pos + 1) {
				cur_pos--;
				if (cur_pos == -1) break;
				part_change[cur_pos] = part_change[cur_pos] + 1;
				for (j = cur_pos + 1; j < num_mismatch; j++)
					part_change[j] = part_change[cur_pos] + j - cur_pos;
			}
			if (cur_pos == -1) break;
		}

		__ASSERT(num_lists == (int)mask_sets.size(), "internal error: wrong num_lists.\n");

		for (i = 0; i < num_lists; i++) {
			int64 mask = mask_sets[i][0];
			__ASSERT((mask & ((1 << num_insdel) - 1)) == 0, "internal error: wrong mask.\n");
			_insdel_masks(mask_sets[i], mask, num_insdel);
		}

		int num_mask_sets = 0;
		int max_num_mask_sets = 0;
		for (i = 0; i < num_lists; i++) {
			num_mask_sets += (int)mask_sets[i].size();
			max_num_mask_sets = max(max_num_mask_sets, (int)mask_sets[i].size());
		}
		if (verbose) cout << "total " << num_mask_sets << " copies of probes.\n";

		used_memory += num_lists * 65536 * 4 * 8; //mask0_table
		used_memory += (size_t)num_probes * num_mask_sets * 4; //index_list
		if (store_key) used_memory += (size_t)num_probes * num_mask_sets * 8; //key_list
		used_memory += (size_t)num_probes * max_num_mask_sets * 16; 
		if (fast_index) used_memory += (size_t)((size_t)num_probes * num_mask_sets * 4 / fast_index_fraction); //fast_indexes
		if (verbose) cout << "estimated memory usage: " << (int)(used_memory >> 20) << "MB.\n";
		if (used_memory > available_memory) {
			panic("not enough memory.\n");
		}

		if (verbose) cout << "Building " << num_lists << " probe lists";
		probe_lists.clear();
		for (i = 0; i < num_lists; i++) {
			if (verbose) cout << "." << flush;
			probe_list *temp_probe_list = new probe_list();
			temp_probe_list->pm = this;
			temp_probe_list->set_mask(mask_sets[i]);
			probe_lists.push_back(temp_probe_list);
		}

		if (verbose) cout << (int)probe_lists.size() << " probe lists created." << endl;
	}

//test speed
	clock_t clock1 = clock();
	i = 0;
	while(true) {
		i++;
		int64 key = (int64)rand() * (int64)rand() * (int64)rand() * (int64)rand();
		correct_key(key);
		vector<uint> indexes;
		search_key(key, indexes);
		clock_t clock2 = clock();
		if (difftime_clock(clock2, clock1) >= 1) break;
	}
	clock_t clock2 = clock();
	if (verbose) cout << "estimated search speed is " << ((search_both_strands && !duplicate_probes)?(i/difftime_clock(clock2, clock1)/2.0):(i/difftime_clock(clock2, clock1))) << " bps/sec.\n";

	return true;
}

inline bool probe_match::initialize(const string &fa_file_name){
	if (!file_exists(fa_file_name)) panic("probe file doesn't exist.\n");
	probe_input_file_name = fa_file_name;
	if (is_fasta_file(fa_file_name)) {
		fasta myfasta;
		myfasta.read_from_file(fa_file_name, false);
		return initialize(myfasta.sequences, myfasta.tags);
	} else if (is_dna_file(fa_file_name)) {
		vector<string> probes;
		probes.clear();
		if (read_dna_from_file(fa_file_name, probes)) {
			return initialize(probes);
		} else {
			panic("probe file is in bad format.\n");
		}
	} else {
		panic("probe file is in bad format.\n");
	}
	return false;
}

inline bool probe_match::initialize(vector<string> &probe_sequences, vector<string> &probe_tags){
	ofstream file2, file3;
	if (output_filtered_probes) {
		file2.open((probe_input_file_name+".filtered").c_str());
		file3.open((probe_input_file_name+".removed").c_str());
	}
	set<string> filter_selected_probe_ids;
	filter_selected_probe_ids.clear();
	if (filter_selected_probes) {
		if (verbose) cout << "reading probe filtering list.\n";
		vector<string> temp_probes;
		if (!read_dna_from_file(filter_selected_probes_filename, temp_probes)) panic("probe filtering file is in bad format.\n");
		for (int i = 0; i < (int)temp_probes.size(); i++) filter_selected_probe_ids.insert(temp_probes[i]);
	}

	if (verbose) cout << "analysing probes...";
	__ASSERT(probe_sequences.size() > 0 && probe_sequences.size() == probe_tags.size(), "internal error, probe_sequences has bad size.\n");
	int min_len = 1000000000, max_len = -1;
	for (int i = 0; i < (int)probe_sequences.size(); i++) {
		probe_sequences[i] = trim_space(toupper(probe_sequences[i]));
		probe_tags[i] = trim_space(probe_tags[i]);
		min_len = min(min_len, (int)probe_sequences[i].length());
		max_len = max(max_len, (int)probe_sequences[i].length());
	}
	probe_len = min_len;
	if (probe_len < 2) probe_len = 2;
	if (probe_len > MAX_PROBE_LEN - num_insdel - 1) probe_len = MAX_PROBE_LEN - num_insdel - 1;
	if (cut_start > 0) probe_len = min(cut_end - cut_start + 1, probe_len);
	if (verbose) cout << "totally " << (int)probe_sequences.size() << " probes, minimum length = " << min_len << ", maximum length = " << max_len << ", set internal key length = " << probe_len << endl;

	if (verbose) cout << "importing probes...";
	int num_probe_skipped = 0;
	for (int i = 0; i < (int)probe_sequences.size(); i++) {
		if (verbose && i % 200000 == 0) cout << "." << flush;
		string probe_id = trim_space(probe_tags[i]);
		probe_info temp_probe_info;
		temp_probe_info.valid = true;
		temp_probe_info.rematch = false;
		temp_probe_info.key = 0;
		temp_probe_info.id = probe_id;
		temp_probe_info.repeat = -1;
		string probe = trim_space(toupper(probe_sequences[i]));
		if (cut_start > 0) {
			probe = probe.substr(cut_start - 1, cut_end - cut_start + 1);
		}

		if (probe.length() == 0) {
			if (verbose) cout << "bad probe length" << endl;
			break;
		}
		temp_probe_info.seq = probe;
		if (filter_selected_probes) {
			if (filter_selected_probe_ids.count(probe_id)) {
				if (output_filtered_probes) file3 << ">" << probe_id << endl << probe << endl;
				temp_probe_info.valid = false;
			}
		}
		if ((int)probe.length() != probe_len) {
			if ((int)probe.length() > probe_len) {
				probe = probe.substr(0, probe_len);
				temp_probe_info.rematch = true;
			} else if ((int)probe.length() < probe_len && (match_shorter_probes)) {
				temp_probe_info.rematch = true;
			} else {
				if (verbose && num_probe_skipped < 10) cout << "Error reading probe " << probe_id << ":" << temp_probe_info.seq << ". probably too short, skipped." << endl;
				num_probe_skipped++;
				if (output_filtered_probes) file3 << ">" << probe_id << endl << temp_probe_info.seq << endl;
				temp_probe_info.valid = false;
			}
		}
		probe = toupper(probe);
		if (filter_low_quality_probes) {
			if (probe.find("AAAAA") != string::npos || probe.find("CCCCC") != string::npos || probe.find("TTTTT") != string::npos || probe.find("GGGGG") != string::npos) {
				if (output_filtered_probes) file3 << ">" << probe_id << endl << temp_probe_info.seq << endl;
				temp_probe_info.valid = false;
			}
		}
		if (!temp_probe_info.valid) {
			probes.push_back(temp_probe_info);
			continue;
		}
		if (match_shorter_probes && (int)probe.length() < probe_len) {
			fix_short_probe(probe, num_mismatch);
			temp_probe_info.rematch = true;
		}
		if (match_N) {
			if (fix_N_probe(probe, num_mismatch)) {
				temp_probe_info.rematch = true;
			}
		}
		__ASSERT((int)probe.length() == probe_len, "internal error, bad probe.\n");
		temp_probe_info.key = tonumber(probe);
		if (temp_probe_info.key != -1) {
			temp_probe_info.valid = true;
			if (output_filtered_probes) file2 << ">" << probe_id << endl << temp_probe_info.seq << endl;
		} else {
			if (verbose && num_probe_skipped < 10) cout << "bad probe charactor in probe " << probe_id << ":" << temp_probe_info.seq << ". skipped\n";
			num_probe_skipped++;
			if (output_filtered_probes) file3 << ">" << probe_id << endl << temp_probe_info.seq << endl;
			temp_probe_info.valid = false;
		}
		probes.push_back(temp_probe_info);
	}
	if (verbose && num_probe_skipped > 0) cout << num_probe_skipped << " probes skipped.\n";
	num_probes = (int)probes.size();
	if (output_filtered_probes) {
		file2.close();
		file3.close();
	}

	if (cluster_identical_probes) {
		if (verbose) cout << "detecting repeated probes...";
		int num_repeat_probes = 0;
		vector<pair<string, int> > temp_probes;
		temp_probes.resize(num_probes);
		for (int i = 0; i < num_probes; i++) temp_probes[i] = pair<string, int>(probes[i].seq, i);
		sort(temp_probes.begin(), temp_probes.end());
		int first_probe_id = -1;
		for (int i = 0; i < num_probes; i++) {
			if (first_probe_id < 0 || temp_probes[i].first != probes[first_probe_id].seq) {
				first_probe_id = temp_probes[i].second;
				probes[temp_probes[i].second].repeat = -1;
			} else {
				probes[temp_probes[i].second].repeat = first_probe_id;
				probes[temp_probes[i].second].valid = false;
				num_repeat_probes++;
			}
		}
		if (verbose) cout << "found " << num_repeat_probes << " repeated probes.\n";
	}

	return initialize();
}

inline bool probe_match::search_key(const int64 key, vector<uint> &indexes){
	vector<uint> temp1, temp2;
	temp2.clear();
	int i, j;

	if (use_hash) {
		int ll;
		if (hash_detect_probe(key, temp1) != (uint)(-1)) {
			for (ll = 0; ll < (int)temp1.size(); ll++) {
				temp2.push_back(decode_probe_index(temp1[ll]));
			}
		}

		if (num_mismatch >= 1 && !hash_1bp_mismatch || num_mismatch >= 2) {
			int num_change = hash_1bp_mismatch ? num_mismatch - 1 : num_mismatch;
			int pos_change[max_mismatch - 1];
			int char_change[max_mismatch - 1];
			int jj;
			for (jj = 0; jj < num_change; jj++) {
				pos_change[jj] = 0;
				char_change[jj] = 0;
			}
			while (true) {
				int64 temp = key;
				for (jj = 0; jj < num_change; jj++)
					temp = change_probe(temp, pos_change[jj], char_change[jj]);
				if (hash_detect_probe(temp, temp1) != (uint)(-1)) {
					for (ll = 0; ll < (int)temp1.size(); ll++) {
						temp2.push_back(decode_probe_index(temp1[ll]));
					}
				}
				int cur_pos = num_change - 1;				
				char_change[cur_pos] = char_change[cur_pos] + 1;
				while (char_change[cur_pos] == 4) {
					char_change[cur_pos] = 0;
					cur_pos--;
					if (cur_pos == -1) break;
					char_change[cur_pos] = char_change[cur_pos] + 1;
				}
				if (cur_pos == -1) {
					cur_pos = num_change - 1;				
					pos_change[cur_pos] = pos_change[cur_pos] + 1;
					while (pos_change[cur_pos] == probe_len) {
						cur_pos--;
						if (cur_pos == -1) break;
						pos_change[cur_pos] = pos_change[cur_pos] + 1;
						for (jj = cur_pos + 1; jj < num_change; jj++)
							pos_change[jj] = pos_change[cur_pos];
					}
				}
				if (cur_pos == -1) break;
			}
		}
	} else {
		for (i = 0; i < (int)probe_lists.size(); i++) {
			temp1.clear();
			if (probe_lists[i]->search_probe(key, temp1)){
				for (j = 0; j < (int)temp1.size(); j++) {
					temp2.push_back(temp1[j]);
				}
			}
		}
	}

	indexes.clear();
	sort(temp2.begin(), temp2.end());
	for (j = 0; j < (int)temp2.size(); j++) {
		if (j == 0 || (temp2[j] != temp2[j - 1])) {
			indexes.push_back(temp2[j]);
		}
	}
	return (indexes.size() > 0);
}

inline void probe_match::search_results_prehandler(vector<search_result_struct> &results, search_result_handler handler) {
	int i;
	vector<search_result_struct> temp_results;
	temp_results.clear();

	for (i = 0; i < (int)results.size(); i++) {
		if (!results[i].reported) {
			temp_results.push_back(results[i]);
		}
	}

	handler(temp_results);
}

inline bool probe_match::search_file(const string &fa_file_name,  search_result_handler handler){
	if (!file_exists(fa_file_name)) panic("transcripts file doesn't exist.\n");
	int i, j, k;
	clock_t clock1 = clock();
	if (verbose) cout << "detecting probes" << endl;
	ifstream file2(fa_file_name.c_str());

	transcripts.clear();
	int current_transcript = 0;
	string transcript_id = "";
	string next_trans_id = "";
	int64 basepair_processed = 0;
	string line = "";
	int64 key = 0;
	int pos = zero_indexed?0:1;
	int last_N = pos-1;
	int current_pos = 0;
	int key_len = probe_len + num_insdel;
	int len_buf = 512;
	vector<search_result_struct> search_results;
	search_results.clear();
	int num_bad_chr_skipped = 0;
 	while (true) {
		while (current_pos > (int)line.length() - len_buf) {
			if (next_trans_id != "") break;
			string templine;
			if (getline(file2, templine)) {
				if (templine[0] == '>') {
					next_trans_id = templine.substr(1, string::npos);
					break;
				} else {
					line = line + templine;
				}
			} else break;
		}
		if (current_pos >= (int)line.length() + num_insdel * 2) {
			if (search_results.size() > 0) search_results_prehandler(search_results, handler);
			search_results.clear();

			if (next_trans_id == "") break;

			transcript_id = next_trans_id;
			transcripts.push_back(transcript_id);
			next_trans_id = "";
			current_transcript++;
//			if (verbose && current_transcript > 1) cout << endl;
//			if (verbose) cout << "processing transcript " << transcript_id;
			pos = zero_indexed?0:1;
			last_N = pos-1;
			current_pos = 0;
			line = "";
			key = 0;
			continue;
		}
		char ch;
		if (current_pos < (int)line.length()) ch = line[current_pos];
		else ch = random_NT();

		if (ch >= 'a' && ch <= 'z') {
			if (search_repeats) ch = (char)toupper(ch);
			else ch = 'N';
		}
		if (ch == 'N') {
			last_N = pos;
			ch = 'A';
		}
		for (i = 0; i < 4; i++) {
			if (ch == DNA_C[i]) break;
		}

		if (i >= 4) {
			if (verbose && num_bad_chr_skipped < 10) cout << "Bad charactor " << ch << " found when processing transcript " << transcript_id << ". Skipped." << endl;
			num_bad_chr_skipped++;
			last_N = pos;
			ch = 'A';
			i = 0;
		}
		key <<= 2;
		correct_key(key, key_len);
		key |= i;
		current_pos++;
		while (current_pos > len_buf + 1000) {
			current_pos -= 1000;
			line = line.substr(1000, string::npos);
		}
		pos++;
		if (verbose && pos % (1000000) == 0) cout << "."<< flush;
		basepair_processed++;
		if (pos - last_N <= key_len) continue;
		vector<uint> indexes;
		indexes.clear();
		vector<search_result_struct> results;
		results.clear();
		if (search_key(key, indexes)) {
			for (j = 0; j < (int)indexes.size(); j++) {
				search_result_struct temp_result;
				temp_result.trans_key_len = key_len;
				temp_result.probe_index = indexes[j];
				if (probes[indexes[j]].rematch || current_pos > (int)line.length()) {
					string probe_seq = probes[indexes[j]].seq;
					int temp_len = (int)probe_seq.length();
					int start = current_pos - key_len;
					int end = min((int)line.length(), current_pos - key_len + temp_len + num_insdel);
					if (duplicate_probes && (int)indexes[j] >= num_probes / 2) {
						probe_seq = get_reverse2(probe_seq);
					}
					__ASSERT(start >= 0 && start < end && end <= (int)line.length(), "internal error: bad trans_seq\n");
					temp_result.trans_seq = toupper(line.substr(start, end - start));
					temp_result.trans_key_len = end - start;
					temp_result.num_mismatch = mismatch(temp_result.trans_seq, probe_seq, num_mismatch, num_insdel, temp_result.trans_key_len);
					temp_result.trans_seq = temp_result.trans_seq.substr(0, temp_result.trans_key_len);
					temp_result.trans_key = key;
				} else {
					temp_result.num_mismatch = mismatch(key, probes[indexes[j]].key, num_mismatch, num_insdel, temp_result.trans_key_len);
					__ASSERT(temp_result.num_mismatch <= num_mismatch, "internal errro: num_mismatch too big\n");
					temp_result.trans_key = (key >> (2 * (key_len - temp_result.trans_key_len)));
					temp_result.trans_seq = tostring(temp_result.trans_key, temp_result.trans_key_len);
				}
				if (temp_result.num_mismatch > num_mismatch) continue;
				temp_result.trans_id = current_transcript - 1;
				temp_result.trans_coord = pos - key_len;
				temp_result.reverse_strand = false;
				if (duplicate_probes && temp_result.probe_index >= num_probes / 2) {
					temp_result.reverse_strand = true;
					temp_result.probe_index -= num_probes / 2;
				}

				temp_result.reported = false;
				results.push_back(temp_result);
			}
		}
		indexes.clear();
		if (search_both_strands && (!duplicate_probes) && search_key(get_reverse(key, key_len), indexes)) {
			for (j = 0; j < (int)indexes.size(); j++) {
				search_result_struct temp_result;
				temp_result.trans_key_len = key_len;
				temp_result.probe_index = indexes[j];
				if (probes[indexes[j]].rematch || current_pos > (int)line.length()) {
					string probe_seq = probes[indexes[j]].seq;
					int temp_len = (int)probe_seq.length();
					int start = max(0, current_pos - temp_len - num_insdel);
					int end = min((int)line.length(), current_pos);
					probe_seq = get_reverse2(probe_seq);
					__ASSERT(start >= 0 && start < end && end <= (int)line.length(), "internal error: bad trans_seq\n");
					temp_result.trans_seq = toupper(line.substr(start, end - start));
					temp_result.trans_key_len = end - start;
					temp_result.num_mismatch = mismatch(temp_result.trans_seq, probe_seq, num_mismatch, num_insdel, temp_result.trans_key_len);
					temp_result.trans_seq = temp_result.trans_seq.substr(0, temp_result.trans_key_len);
					temp_result.trans_key = key;
				} else {
					temp_result.num_mismatch = mismatch(get_reverse(key, key_len), probes[indexes[j]].key, num_mismatch, num_insdel, temp_result.trans_key_len);
					__ASSERT(temp_result.num_mismatch <= num_mismatch, "internal errro: num_mismatch too big\n");
					temp_result.trans_key = (key >> (2 * (key_len - temp_result.trans_key_len)));
					temp_result.trans_seq = tostring(temp_result.trans_key, temp_result.trans_key_len);
				}
				if (temp_result.num_mismatch > num_mismatch) continue;
				temp_result.trans_id = current_transcript - 1;
				temp_result.trans_coord = pos - key_len;
				if (probes[indexes[j]].rematch) temp_result.trans_coord = pos - (int)probes[indexes[j]].seq.length() - num_insdel;
				temp_result.reverse_strand = true;

				temp_result.reported = false;
				results.push_back(temp_result);
			}
		}

		if (filter_results) {
			vector<search_result_struct> temp_results;
			temp_results.clear();
			for (j = 0; j < (int)search_results.size(); j++) {
				bool remain = true;
				for (k = 0; k < (int)results.size(); k++) {
					if (filter_results && search_results[j].probe_index == results[k].probe_index && search_results[j].num_mismatch > results[k].num_mismatch) {
						remain = false;
						break;
					}
				}
				if (remain && !search_results[j].reported) temp_results.push_back(search_results[j]);
			}
			if (temp_results.size() > 0) search_results_prehandler(temp_results, handler);
			temp_results.clear();
			for (j = 0; j < (int)results.size(); j++) {
				for (k = 0; k < (int)search_results.size(); k++) {
					if (filter_results && search_results[k].probe_index == results[j].probe_index && search_results[k].num_mismatch <= results[j].num_mismatch) {
						results[j].reported = true;
						break;
					}
				}
			}
			search_results.clear();
			search_results = results;
		} else {
			handler(results);
		}
		results.clear();
	}

	file2.close();
	clock_t clock2 = clock();
	if (verbose && num_bad_chr_skipped > 0) cout << num_bad_chr_skipped << " bad characters skipped.\n";
	if (verbose) cout << endl;
#ifdef __DEBUG
	if (verbose && (!use_hash)) cout << "average search steps:" << (double)num_total_search_step / num_search << " maximum search steps:" << num_max_search_step << endl;
#endif
	if (verbose) cout << basepair_processed << " base pairs processed. average search speed: " << basepair_processed / difftime_clock(clock2, clock1) << " bps/sec.\n";
//	if (verbose) cout << "memory used: " << (int)(used_memory >> 20) << "MB.\n";
	return true;
}

#endif//PROBE_MATCH_H
