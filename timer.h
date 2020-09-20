#ifndef TIMER_H
#define TIMER_H

#include "stl.h"

inline double difftime_clock(clock_t &clock2, clock_t &clock1) {
	return double(clock2 - clock1)/CLOCKS_PER_SEC;
}

class timer{
public:
	bool auto_summary;
	bool auto_click;
	bool accumulate_time;
	bool auto_report;
	bool use_clock;

	vector<time_t> times;
	vector<clock_t> clocks;

	timer(){
		use_clock = true;
		auto_summary = true;
		auto_click = true;
		accumulate_time = false;
		auto_report = false;
		reset();
	}
	~timer(){
		if (auto_click) click();
		if (auto_summary) print_usage();
	}
	void reset(){
		times.clear();
		clocks.clear();
		click();
	}
	void click(){
		if (use_clock) clocks.push_back(clock());
		else times.push_back(time((time_t*) NULL));
		if (auto_report) print_usage();
	}
	void print_usage(){
		int i;
		if (use_clock) {
			for (i = 1; i < (int)clocks.size(); i++) {
				if (accumulate_time)
					cout << "time used: " << difftime_clock(clocks[i], clocks[0]) << " seconds\n";
				else 
					cout << "time used: " << difftime_clock(clocks[i], clocks[i-1]) << " seconds\n";
			}
		} else {
			for (i = 1; i < (int)times.size(); i++) {
				if (accumulate_time)
					cout << "time used: " << difftime(times[i], times[0]) << " seconds\n";
				else 
					cout << "time used: " << difftime(times[i], times[i-1]) << " seconds\n";
			}
		}
	}
};

#endif //TIMER_H
