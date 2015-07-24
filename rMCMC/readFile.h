/*
 * readFile.h
 *
 *  Created on: Oct 8, 2014
 *      Author: claudebherer
 */

#ifndef READFILE_H_
#define READFILE_H_



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include "event.h"

using namespace std;

void readFile(string &filename, vector<event> &events);



#endif /* READFILE_H_ */
