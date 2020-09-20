#ifndef MEMORY_UTILS_H
#define MEMORY_UTILS_H

#include "type.h"
#include "stl.h"
#include "timer.h"

inline size_t get_available_memory_size(bool timed_operation = true){
	vector<char*> ptrs;
	ptrs.clear();
	double time_interval = 0;
	clock_t clock0 = clock();
	while(true) {
		clock_t clock1 = clock();
		char* ptr = (char*)malloc(1 << 26);
		if (ptr == NULL) break;
		memset(ptr, -1, 1 << 26);
		ptrs.push_back(ptr);
		clock_t clock2 = clock();
		if (ptrs.size() == 1) time_interval = difftime_clock(clock2, clock1);
		if (timed_operation && (difftime_clock(clock2, clock1) >= time_interval * 3 || difftime_clock(clock2, clock0) >= 10)) break;
	}
	size_t result = ptrs.size() << 26;
	for (int i=0; i < (int)ptrs.size(); i++) delete[] ptrs[i];
	ptrs.clear();
	return result;
}

inline bool is32system(){
	return (4 == sizeof(void*));
}

inline bool is_little_endian(){
	uint x = 1;
	return (1 == *((char*)(&x)));
}

#endif //MEMORY_UTILS_H
