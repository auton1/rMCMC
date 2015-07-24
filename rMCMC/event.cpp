/*
 * events.cpp
 *
 *  Created on: Oct 9, 2014
 *      Author: claudebherer
 */

#include "event.h"

void event::set_seed(int seed)
{
	generator.seed(seed);
}

void event::draw_new_current_interval(vector< interval > &intervals)
{
	static vector<double> lambda_subset;	// Keep as static to avoid repeatedly reserving memory
	lambda_subset.resize(N_intervals);
	double sum = 0;
	for (unsigned int ui=0; ui<N_intervals; ui++)
	{
		lambda_subset[ui] = intervals[min_interval_idx+ui].lambdaij;
	}

	discrete_distribution<unsigned int> distribution(lambda_subset.begin(), lambda_subset.end());
	current_interval_idx = min_interval_idx + distribution(generator);
}

void event::take_interval_sample()
{
	interval_samples[current_interval_idx - min_interval_idx]++;
	if (interval_samples[current_interval_idx - min_interval_idx] == 0)
	{
		cerr << "Error: Overflow! interval_samples cannot be of type unsigned short!" << endl;
		exit(1);
	}
}
