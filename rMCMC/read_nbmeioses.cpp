/*
 * read_nbmeioses.cpp
 *
 *  Created on: Nov 20, 2014
 *      Author: claudebherer
 */


#include "read_nbmeioses.h"


// Read nb_meioses file to a vector of events
// this file includes the start and end positions of intervals, and the nb of meioses
// format: start end nbmeioses

void read_nbmeioses(string &M_filename, vector<interval> &intervals)
{
  int pos1, pos2;
  int nbmeioses;
  string line;
  stringstream ss;

  ifstream myfile(M_filename.c_str());
  if (!myfile.good())
  {
	  // Error
	  cerr << "cannot open file" << endl;
	  exit(1);
  }

  // Build the vector of intervals while reading file
  intervals.resize(0);
  interval my_interval;

  while(!myfile.eof())
  {
	  getline(myfile, line);
	  if (line == "")
		  continue;

	  ss.clear(); ss.str()="";
	  ss.str(line);
	  ss >> pos1 >> pos2 >> nbmeioses;

	  my_interval.startpos = pos1;
	  my_interval.stoppos = pos2;
	  my_interval.M = nbmeioses;

	  intervals.push_back(my_interval);

  }
  myfile.close();
}; // fin read_nbmeioses
