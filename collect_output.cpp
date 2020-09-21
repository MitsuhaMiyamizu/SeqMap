#include "stl.h"
#include "string_operation.h"
#include "file_operation.h"

class mini_search_result_struct{
public:
	string trans_id;
	int trans_coord;
	int num_mismatch;
};

bool operator< (const mini_search_result_struct &result1, const mini_search_result_struct &result2) {
	return result1.num_mismatch < result2.num_mismatch;
}

class probe_info{
public:
	vector<int> counts;
	vector<mini_search_result_struct> results;
};

map<string, probe_info> results;

int main(int argc, char* argv[]){
	int i;
	if (argc < 2) {
		cout << "usage: collect_output <files_list or filenames> \n";
		return 1;
	}
	ifstream list_file;
	if (argc == 2) {
		string list_filename = argv[1];
		list_file.open(list_filename.c_str());
	}
	results.clear();
	int file_i = 0;
	int num_mismatch = -1;
	int top_results = -1;
	bool mismatch_pos = false;
	while(true) {
		file_i++;
		string line;
		if (argc == 2) {
			getline(list_file, line);
			if (line == "") break;
		} else {
			if (file_i == argc) break;
			line = argv[file_i];
		}
		if (!file_exists(line)) panic(string("file does not exist:") + line);
		cout << "processing file: " << line << endl;
		ifstream file(line.c_str());
		bool first_line = true;
		while (true) {
			string read_line;
			getline(file, read_line);
			if (read_line == "") break;
			if (first_line) {
				size_t start = 0;
				int n_mismatch = -1;
				while ((start = read_line.find("#mismatch=", start)) != string::npos) {
					n_mismatch++;
					start += 1;
				}
				if (num_mismatch == -1) num_mismatch = n_mismatch;
				else __ASSERT(num_mismatch == n_mismatch, "inconsistent num_mismatch.\n");

				start = 0;
				int n_top_results = 0;
				while ((start = read_line.find("coord", start)) != string::npos) {
					n_top_results++;
					start += 1;
				}
				if (top_results == -1) top_results = n_top_results;
				else __ASSERT(top_results == n_top_results, "inconsistent top_results.\n");

				mismatch_pos = (read_line.find("mismatch_pos") != string::npos);

				first_line = false;
				continue;
			}
			char probe_id[1024];
			probe_info info;
			int temp;
			vector<string> tokens = string_tokenize(read_line, " \t");
			if (1 != sscanf(tokens[0].c_str(), " %s", probe_id)) panic("wrong format.\n");
			for (i = 0; i <= num_mismatch; i++) {
				if (1 != sscanf(tokens[i+1].c_str(), " %d", &temp)) panic("wrong format.\n");
				info.counts.push_back(temp);
			}
			
			for (i = 0; i <= top_results; i++) {
				int index = num_mismatch + 2 + i * 3;
				if (mismatch_pos) index++;
				if (index >= (int)tokens.size()) break;
				mini_search_result_struct result;
				char trans_id[1024];
				if (1 != sscanf(tokens[index].c_str(), " %s", trans_id)) panic("wrong format.\n");
				result.trans_id = trans_id;
				if (1 != sscanf(tokens[index+1].c_str(), " %d", &temp)) panic("wrong format.\n");
				result.trans_coord = temp;
				if (1 != sscanf(tokens[index+2].c_str(), " %d", &temp)) panic("wrong format.\n");
				result.num_mismatch = temp;
				info.results.push_back(result);
			}
			if (results.count(probe_id) == 0) {
				results[probe_id] = info;
			} else {
				probe_info old_info = results[probe_id];
				for (i = 0; i <= num_mismatch; i++) old_info.counts[i] += info.counts[i];
				for (i = 0; i < (int)info.results.size(); i++) old_info.results.push_back(info.results[i]);
				results[probe_id] = old_info;
			}
		}
		file.close();
	}
	if (argc == 2) {
		list_file.close();
	}

	cout << "writing results.\n";
	ofstream output_file("all.output");
	output_file << "probe_id";
	for (i = 0; i <= num_mismatch; i++) output_file << "\t#mismatch=" << i;
	for (i = 0; i < top_results; i++) output_file << "\ttrans_id\tcoord\t#mismatch";
	output_file << endl;
	map<string, probe_info>::iterator map_iterator;
	for (map_iterator = results.begin(); map_iterator != results.end(); map_iterator++) {
		output_file << map_iterator->first;
		for (i = 0; i <= num_mismatch; i++) output_file << "\t" << (map_iterator->second.counts)[i];
		sort(map_iterator->second.results.begin(), map_iterator->second.results.end());
		for (i = 0; i < min(top_results, (int)map_iterator->second.results.size()); i++) output_file << "\t" << (map_iterator->second.results)[i].trans_id << "\t" << (map_iterator->second.results)[i].trans_coord << "\t" << (map_iterator->second.results)[i].num_mismatch;
		output_file << endl;
	}
	output_file.close();
	return 0;
}
