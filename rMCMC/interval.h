/*
 * interval.h
 *
 *  Created on: Oct 9, 2014
 *      Author: claudebherer
 */

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

class interval
{
public:
	interval(){};	// constructor
	~interval(){}; // destructor

	int startpos, stoppos;
	double r0;
	
	int M; // nb of meioses in this interval

	double lambda0;
	double alpha0, beta0;
	static double prior_recomb_rate;
	static double prior_variance_rate;	// Variance per megabase

	double lambdaij;	// Recombination fraction
	double alphaij, betaij;

	double MLE;

	static int CI95_tresh;
	vector<double> lower95CI;
	vector<double> upper95CI;

	static int CI99_tresh;
	vector<double> lower99CI;
	vector<double> upper99CI;

	void set_gamma_priors();
	void draw_new_lamba();
	void set_alpha_and_beta(double alpha, double beta);
	void set_nb_events_currently(int nb_events_curr);
	void set_nbMeioses(int m);
	void set_MLE(double mle);
	void set_CI(double rate_sample);
	static void set_seed(int seed);

private:
	static std::default_random_engine generator;
};




#endif /* INTERVAL_H_ */
