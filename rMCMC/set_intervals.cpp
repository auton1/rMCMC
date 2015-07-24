/*
 * set_intervals.cpp
 *
 *  Created on: Oct 10, 2014
 *      Author: claudebherer
 */

#include "set_intervals.h"

void set_intervals(vector<event> &events, vector<interval> &intervals)
{
	int nevents;
	nevents = events.size();

 	// For each event, find min_interval_idx and max_interval.idx
 	for (int i=0; i<events.size(); i++)
 	{
 		for (int j=0; j<intervals.size(); j++)
 		{
 			if (events[i].min_interval_pos == intervals[j].startpos)
 				events[i].min_interval_idx = j;
 			if (events[i].max_interval_pos == intervals[j].stoppos)
 				events[i].max_interval_idx = j;
 		}
		//cerr << events[i].min_interval_pos << "\t" << events[i].max_interval_pos << endl;

		events[i].N_intervals = events[i].max_interval_idx - events[i].min_interval_idx + 1;
		events[i].interval_samples.resize(events[i].N_intervals, 0);
 	}


 	/*
	// PRINT the events that overlap each interval i
	cerr << "intervals:\t size:" << intervals.size() << endl;
	for (int i=0; i<intervals.size(); i++){
		cerr << "interval " << i << " :" << endl;
		cerr << intervals[i].startpos << "\t" << intervals[i].stoppos << endl;
		//cout << "d_bp: " << intervals[i].stoppos-intervals[i].startpos << endl;

		vector<int> event_ids;
		event_ids = intervals[i].events_ID;
		cerr << event_ids.size() << " overlapping events : \n";
		for (int j=0; j<event_ids.size(); j++)
		{
			cerr << events[event_ids[j]].min_interval_pos << "\t" << events[event_ids[j]].max_interval_pos <<  endl;
		}
	}
	*/


}

