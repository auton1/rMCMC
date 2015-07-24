/*
 * readFile.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: claudebherer
 */

#include "readFile.h"

// Read event file to a vector of events

void readFile(string &filename, vector<event> &my_vector)
{
  int pos1, pos2;
  string line;
  stringstream ss;

  ifstream myfile(filename.c_str());
  if (!myfile.good())
  {
	  // Error
	  cerr << "cannot open file" << endl;
	  exit(1);
  }

  my_vector.resize(0);
  event my_event;
  while(!myfile.eof())
  {
	  getline(myfile, line);
	  if (line == "")
		  continue;

	  ss.clear(); ss.str()="";
	  ss.str(line);
	  ss >> pos1 >> pos2;

	  my_event.min_interval_pos = pos1;
	  my_event.max_interval_pos = pos2;
	  my_vector.push_back(my_event);
  }
  myfile.close();
}; // fin readFile


/* Print the content of my vector
for (int i=0;i<row;i++) {
  cout << i << ": " ;
  for (int j=0;j<2;j++) {
    cout << "\t" << my_vector[i][j];
  }
  cout << endl;
}*/
