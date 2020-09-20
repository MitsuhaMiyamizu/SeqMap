#include "probe_match.h"
#include "string_operation.h"
#include "timer.h"

probe_match myprobe_match;
ofstream file;
bool statistics_only = true;
bool output_probe_without_match = true;
bool output_alignment = false;
int num_top_matches = 2;
int num_probes;
bool do_eland = true;
int eland_version = 2;

class mini_search_result_struct{
public:
	int trans_id;
	int trans_coord;
	bool reverse_strand;
	uchar num_mismatch;
};

int* match_counts;
vector<vector<mini_search_result_struct> > top_matches;

void file_write_handler(vector<search_result_struct> &results) {
	int i, j;
	for (i = 0; i < (int)results.size(); i++) {
		__ASSERT(!results[i].reported, "internal error: reported.\n");		
		if (statistics_only) {
			__ASSERT(results[i].probe_index >= 0 && results[i].probe_index < num_probes, "internal error: wrong probe_index.\n");
			__ASSERT(results[i].num_mismatch >= 0 && results[i].num_mismatch <= myprobe_match.num_mismatch, "internal error: wrong num_mismatch.\n");
			match_counts[results[i].probe_index * (myprobe_match.num_mismatch + 1) + results[i].num_mismatch]++;
			j = 0;
			for (j = 0; j < (int)top_matches[results[i].probe_index].size(); j++) {
				if (results[i].num_mismatch < top_matches[results[i].probe_index][j].num_mismatch) {
					break; 
				}
			}
			if (j < num_top_matches) {
				mini_search_result_struct result;
				result.trans_id = results[i].trans_id;
				result.trans_coord = results[i].trans_coord;
				result.reverse_strand = results[i].reverse_strand;
				result.num_mismatch = (uchar)results[i].num_mismatch;
				top_matches[results[i].probe_index].insert(top_matches[results[i].probe_index].begin() + j, result);
			}

			if ((int)top_matches[results[i].probe_index].size() > num_top_matches) {
				top_matches[results[i].probe_index].erase(top_matches[results[i].probe_index].begin() + num_top_matches);
				__ASSERT((int)top_matches[results[i].probe_index].size() == num_top_matches, "internal error: wrong num_top_matches.\n");
			}
		} else {
			if (results[i].reverse_strand) results[i].trans_seq = get_reverse2(results[i].trans_seq);
			string trans_seq = results[i].trans_seq, probe_seq = myprobe_match.probes[results[i].probe_index].seq;
			if (output_alignment){
				local_alignment(trans_seq, probe_seq, myprobe_match.num_insdel);
			}
			file << myprobe_match.transcripts[results[i].trans_id] << "\t" << results[i].trans_coord << "\t" << trans_seq;
			file << "\t" << myprobe_match.probes[results[i].probe_index].id;
			file << "\t" << probe_seq << "\t" << results[i].num_mismatch;
			if (myprobe_match.search_both_strands) {
				if (results[i].reverse_strand) file << "\t-";
				else file << "\t+";
			}
			file << endl;
		}
	}
}

void do_statistics(){
	for (int i = 0; i < num_probes; i++) {
		if (myprobe_match.probes[i].repeat >= 0) {
			int mirror_probe = myprobe_match.probes[i].repeat;
			for (int j = 0; j <= myprobe_match.num_mismatch; j++) {
				match_counts[i * (myprobe_match.num_mismatch + 1) + j] = match_counts[mirror_probe * (myprobe_match.num_mismatch + 1) + j];
			}
			top_matches[i] = top_matches[mirror_probe];
		}
	}

	if (do_eland) {
		for (int i = 0; i < num_probes; i++) {
			if (eland_version > 1) file << myprobe_match.probes[i].id << "\t";
			file << myprobe_match.probes[i].seq;
			int top_target = -1;
			if (top_matches[i].size() > 0) top_target = top_matches[i][0].num_mismatch;
			bool unique_target = false;
			if (top_matches[i].size() == 1 || top_matches[i].size() > 1 && top_matches[i][0].num_mismatch < top_matches[i][1].num_mismatch) unique_target = true;
			int* target_counts = match_counts + i * (myprobe_match.num_mismatch + 1);
			if (eland_version == 1) {
				file << "\t" << top_target + 1 << "\t";
				if (top_target >= 0) file << (int)target_counts[top_target]; else file << "0";
			} else if (eland_version == 2) {
				file << "\t";
				if (top_target >= 0) {
					__ASSERT(target_counts[top_target] > 0, "internal error, target_counts[top_target] == 0");
					if (target_counts[top_target] > 1) file << "R"; else file << "U";
					file << top_target;
				} else file << "NM";
				for (int j = 0; j <= myprobe_match.num_mismatch; j++) file << "\t" << (int)target_counts[j];
			} else if (eland_version == 3){
				file << "\t";
				if (top_target < 0) file << "NM";
				else {
					for (int j = 0; j <= myprobe_match.num_mismatch; j++) {
						if (j > 0) file << ":";
						file << (int)target_counts[j];
					}
				}
			} else panic("bad eland_version");
			if (top_target >= 0) {
				if (unique_target && eland_version < 3) {
					file << "\t" << myprobe_match.transcripts[top_matches[i][0].trans_id];
					if (eland_version == 1) file << ":"; else file << "\t";
					file << top_matches[i][0].trans_coord << "\t" << (top_matches[i][0].reverse_strand?"R":"F");
				} else if (eland_version == 3) {
					file << "\t";
					for (int j = 0; j < (int)top_matches[i].size(); j++) {
						if (top_matches[i][j].num_mismatch > top_target + 1) break;
						if (j > 0) file << ",";
						file << myprobe_match.transcripts[top_matches[i][j].trans_id] << ":" << top_matches[i][j].trans_coord << (top_matches[i][j].reverse_strand?"R":"F") << (int)top_matches[i][j].num_mismatch;
					}
				}
			}
			file << endl;
		}
	} else {
		int i, j;
		file << "probe_id";
		for (j = 0; j <= myprobe_match.num_mismatch; j++) {
			file << "\t#mismatch=" << j;
		}
		for (j = 0; j < num_top_matches; j++) {
			file << "\ttrans_id\tcoord\t#mismatch";
		}
		file << endl;
		for (i = 0; i < num_probes; i++) {
			int sum = 0;
			for (j = 0; j <= myprobe_match.num_mismatch; j++) {
				sum += match_counts[i * (myprobe_match.num_mismatch + 1) + j];
			}
			if (sum == 0 && !output_probe_without_match) continue;
			file << myprobe_match.probes[i].id;
			for (j = 0; j <= myprobe_match.num_mismatch; j++) {
				file << "\t" << match_counts[i * (myprobe_match.num_mismatch + 1) + j];
			}
			for (j = 0; j < (int)top_matches[i].size(); j++) {
				file << "\t" << myprobe_match.transcripts[top_matches[i][j].trans_id] << "\t" << top_matches[i][j].trans_coord << "\t" << (int)(top_matches[i][j].num_mismatch);
			}
			file << endl;
		}
	}
}

int main(int argc, char* argv[]){
	timer mytimer;
	int i, j;
	probe_len = -1;
	
	if (argc > 1 && argc < 5) {
		panic("not enough parameters.\n");
	} else if (argc <= 1) {
		printf("SeqMap: Short Sequence Mapping Tool\n\n");
		printf("Usage: seqmap <number of mismatches> <probe FASTA file name> <transcript FASTA file name> <output file name> [options]\n\n");
		printf("parameters/Options: (* for advanced users only)\n");
		printf("<number of mismatches>                         \tmaximum edit distance allowed\n");
		printf("<probe FASTA file name>                        \tprobe/tag/read sequences\n");
		printf("<transcript FASTA file name>                   \treference sequences\n");
		printf("<output file name>                             \tname of the output file\n");
		printf("[/eland[:style]]                               \toutput in eland style\n");
		printf("[/output_top_matches:num_top_matches]          \toutput the top\n\t\t\t\t\t\t \"num_top_matches\" targets\n");
		printf("[/forward_strand]                              \tsearch forward strand only\n");
		printf("[/allow_insdel:num_insdel]                     \tenable insertion and deletion\n");
		printf("[/cut:start,end]                               \tcut the probes\n");
		printf("[/match_shorter_probes]                        \tmatch shorter probes\n");
		printf("[/skip_N]                                      \tskip probes that have N or .\n");
		printf("[/no_repeats]                                  \tdo not search repeat regions\n");
		printf("[/silent]                                      \tdisable debug information\n");
		printf("[/available_memory:memory_size(in MB)]         \tdisable memory detection and\n\t\t\t\t\t\t set memory size manually\n");
		printf("[/zero_indexed]                                \toutput coordinates in\n\t\t\t\t\t\t 0-indexed manner\n");
		printf("[/output_all_matches]                          \toutput all matches\n");
		printf("[/exact_mismatch]                              \toutput targets with\n\t\t\t\t\t\t exact number of mismatches only\n");
		printf("*[/output_alignment]                            \toutput alignment\n");
		printf("*[/output_statistics]                          \toutput in probe order\n");
		printf("*[/do_not_output_probe_without_match]          \tignore non-matched probes\n");
		printf("*[/no_duplicate_probes]                        \tdo not duplicate probes for\n\t\t\t\t\t\t reverse strand search\n");
		printf("*[/use_hash]                                   \tuse hash if possible\n");
		printf("*[/no_hash_1bp_mismatch]                       \tdo not use hash for 1bp search\n");
		printf("*[/limit_num_part:max_num_part]                \tmaximum number of parts that\n\t\t\t\t\t\t a probe is splitted\n");
		printf("*[/interpolation_search]                       \tuse interpolation search\n\t\t\t\t\t\t ,not binary search\n");
		printf("*[/shift_mask]                                 \tmay reduce memory usage\n\t\t\t\t\t\t while increase running time\n");
		printf("*[/no_store_key]                               \tmay reduce memory usage\n\t\t\t\t\t\t while increase running time\n");
		printf("*[/no_fast_index]                              \tmay reduce memory usage\n\t\t\t\t\t\t while increase running time\n");
		printf("*[/fast_index_fraction:fraction]\n");
		printf("*[/no_filter_results]\n");
		printf("*[/filter_selected_probes:selected_probes_file]\tremove the probes in file\n");
		printf("*[/output_filtered_probes]                     \toutput the filtered probes\n");
		printf("*[/filter_low_quality_probes]                  \tfilter low quality probes\n");
		printf("\nVisit software homepage http://biogibbs.stanford.edu/~jiangh/SeqMap/ for more information. Bugs report to jiangh@stanford.edu.\n");
		return 1;
	}	

	for (i = 5; i < argc; i++) {
		string option = argv[i];
		string command = "";
		string parameter = "";
		if (option[0] != '/') {
			panic(string("bad option ") + option);
		}
		size_t index = option.find (":");
		if (index != string::npos ) {
			command = option.substr(1, index - 1);
			parameter = option.substr(index + 1);
		} else {
			command = option.substr(1);
		}

		if (command == "forward_strand") {
			myprobe_match.search_both_strands = false;
		} else if (command == "silent") {
			myprobe_match.verbose = false;
		} else if (command == "exact_mismatch") {
			myprobe_match.exact_mismatch = true;
		} else if (command == "output_alignment") {
			output_alignment = true;
		} else if (command == "no_repeats") {
			myprobe_match.search_repeats = false;
		} else if (command == "interpolation_search") {
			myprobe_match.interpolation_search = true;
		} else if (command == "shift_mask") {
			myprobe_match.shift_mask = true;
		} else if (command == "no_store_key") {
			myprobe_match.store_key = false;
		} else if (command == "no_fast_index") {
			myprobe_match.fast_index = false;
		} else if (command == "limit_num_part") {
			myprobe_match.limit_num_part = str2int(parameter);
		} else if (command == "allow_insdel") {
			myprobe_match.num_insdel = str2int(parameter);
		} else if (command == "fast_index_fraction") {
			myprobe_match.fast_index_fraction = str2double(parameter);
		} else if (command == "cut") {
			vector<string> tokens = string_tokenize(parameter, ",");
			if (tokens.size() != 2) panic("bad cut parameters.\n");
			myprobe_match.cut_start = str2int(tokens[0]);
			myprobe_match.cut_end = str2int(tokens[1]);
			if (!(myprobe_match.cut_start >= 1 && myprobe_match.cut_end > myprobe_match.cut_start)) panic("bad cut parameters.\n");
		} else if (command == "no_duplicate_probes") {
			myprobe_match.duplicate_probes = false;
		} else if (command == "use_hash") {
			myprobe_match.use_hash = true;
		} else if (command == "no_filter_results") {
			myprobe_match.filter_results = false;
		} else if (command == "no_hash_1bp_mismatch") {
			myprobe_match.hash_1bp_mismatch = false;
		} else if (command == "output_statistics") {
			do_eland = false;
			num_top_matches = 2;
		} else if (command == "do_not_output_probe_without_match") {
			output_probe_without_match = false;
		} else if (command == "filter_selected_probes") {
			myprobe_match.filter_selected_probes = true;
			myprobe_match.filter_selected_probes_filename = parameter;
		} else if (command == "output_filtered_probes") {
			myprobe_match.output_filtered_probes = true;
		} else if (command == "filter_low_quality_probes") {
			myprobe_match.filter_low_quality_probes = true;
		} else if (command == "output_top_matches") {
			num_top_matches = str2int(parameter);
		} else if (command == "available_memory") {
			myprobe_match.available_memory = ((int64)(1<<20)) * str2int(parameter);
		} else if (command == "eland") {
			if (is_int(parameter)) eland_version = str2int(parameter);	
			if (eland_version < 1 || eland_version > 3) eland_version = 3;
			if (eland_version == 3) num_top_matches = 10; else num_top_matches = 2;
		} else if (command == "match_shorter_probes") {
			myprobe_match.match_shorter_probes = true;
		} else if (command == "skip_N") {
			myprobe_match.match_N = false;
		} else if (command == "zero_indexed") {
			myprobe_match.zero_indexed = true;
		} else if (command == "output_all_matches") {
			statistics_only = false;
			myprobe_match.cluster_identical_probes = false;
		} else {
			panic(string("bad command ") + command);
		}
	}

	sscanf(argv[1], "%d", &myprobe_match.num_mismatch);
	if (myprobe_match.num_mismatch < 0) {
		myprobe_match.num_mismatch = - myprobe_match.num_mismatch;
		myprobe_match.exact_mismatch = true;
	}
	if (myprobe_match.num_mismatch < 0 || myprobe_match.num_mismatch > max_mismatch) {
		panic(string("ERROR: number of mismatches must be between 0 and ") + int2str(max_mismatch));
	}
	if (myprobe_match.num_insdel < 0 || myprobe_match.num_insdel > max_insdel) {
		panic(string("ERROR: number of insdels must be between 0 and ") + int2str(max_insdel));
	}
	if (myprobe_match.num_insdel > myprobe_match.num_mismatch) {
		panic("ERROR: number of insdels can not be greater than num of mismatches.\n\n");
	}

	if (myprobe_match.verbose) {
		cout << "Command line:";
		for (i = 0; i < argc; i++)
			cout << argv[i] << " ";
		cout << endl;
	}

	if (statistics_only) {
		myprobe_match.extra_memory_per_probe = (myprobe_match.num_mismatch + 1) * 4; //match_counts
		myprobe_match.extra_memory_per_probe += sizeof(vector<mini_search_result_struct>) + num_top_matches * sizeof(mini_search_result_struct); //top_matches;
	}

	myprobe_match.initialize(argv[2]);
	if (myprobe_match.duplicate_probes) num_probes = myprobe_match.num_probes / 2; 
	else num_probes = myprobe_match.num_probes;
	match_counts = NULL;

	if (statistics_only) {
		match_counts = new int[num_probes * (myprobe_match.num_mismatch + 1)];
		top_matches.resize(num_probes);
		for (i = 0; i < num_probes; i++){
			for (j = 0; j <= myprobe_match.num_mismatch; j++) match_counts[i * (myprobe_match.num_mismatch + 1) + j] = 0;
			top_matches[i].clear();
		}
	}

	file.open(argv[4]);
	if (!statistics_only) {
		file << "trans_id\ttrans_coord\ttarget_seq\tprobe_id\tprobe_seq\tnum_mismatch";
		if (myprobe_match.search_both_strands) file << "\tstrand";
		file << endl;
	}
	myprobe_match.search_file(argv[3], file_write_handler);
	if (statistics_only) do_statistics();
	file.close();
	if (match_counts != NULL) delete[] match_counts;
}
