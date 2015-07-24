/*
 * event.h
 *
 *  Created on: Oct 9, 2014
 *      Author: claudebherer
 */

#ifndef EVENT_H_
#define EVENT_H_

#include <random>
#include <vector>
#include <iostream>

#include "interval.h"

using namespace std;

class event
{
public:
	unsigned int min_interval_pos;
	unsigned int max_interval_pos;

	unsigned int min_interval_idx;
	unsigned int max_interval_idx;
	
	unsigned int N_intervals;

	unsigned int current_interval_idx;
	
	vector<unsigned short> interval_samples;	// Can handle up to 65,535 counts.

	event(){};	 // constructor
	~event(){};  // destructor

	void draw_new_current_interval(vector< interval > &intervals);
	static void set_seed(int seed);
	void take_interval_sample();

private:
	static std::default_random_engine generator;
};



#endif /* EVENT_H_ */
