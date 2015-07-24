/*
 * interval.cpp
 */

#include "interval.h"

void interval::set_gamma_priors()
{
	double d_Mb = (stoppos - startpos) / 1000000.0; // d_Mb : distance in Mb
	double d_cM = d_Mb * prior_recomb_rate; // d_cM : distance in cM given a prior recombination rate

	// prior variance in the nb of events per interval
	/*
	static double sd0 = sqrt(prior_variance_rate); // sdev of rec. rates in cM/Mb
	double sd0_d = (sd0 * d_Mb);  	// sdev of genetic distance in centiMorgan
	double var0 = sd0_d*sd0_d; 	//prior variance in the nb of events in a seq interval of size d_Mb
	*/
	// Based on HapMap, we expect the variance to increase proportional to (window size)^1.5.
	double var0 = prior_variance_rate * pow(d_Mb, 1.5);
	double sd0_d = sqrt(var0);
	
	// Convert to recombination fraction
	r0 = 0.5*(exp(0.04*d_cM) - 1.0) / (exp(-0.04*d_cM)+1.0);
	double r0_sd = 0.5*(exp(0.04*sd0_d) - 1.0) / (exp(-0.04*sd0_d)+1.0);

	double c = r0_sd*r0_sd;	// Required variance for this interval
	alpha0 = r0*r0 / c;
	beta0 = r0 / c;
	
	// Mean = alpha / beta = (r0*r0/c) / (r0 / c) = r0;
	// Variance = (alpha / beta) / beta = r0 / (r0 / c) = c

	gamma_distribution<double> distribution(alpha0, 1.0/beta0);
	lambda0 = distribution(generator);
	if (lambda0==0)
		lambda0 = numeric_limits<double>::epsilon();
}

void interval::draw_new_lamba()
{
	gamma_distribution<double> distribution(alphaij, 1.0/betaij);
	lambdaij = distribution(generator);
	if (lambdaij == 0)
		lambdaij = numeric_limits<double>::epsilon();
}

void interval::set_alpha_and_beta(double alpha, double beta)
{
	alphaij = alpha;
	betaij = beta;
}

void interval::set_seed(int seed)
{
	generator.seed(seed);
}

void interval::set_nbMeioses(int m)
{
	M = m;
}

void interval::set_MLE(double mle)
{
	MLE = mle;
}

void interval::set_CI(double rate_sample)
{
	lower95CI.push_back(rate_sample);
	sort(lower95CI.begin(), lower95CI.end());
	if(lower95CI.size()>CI95_tresh)
		lower95CI.resize(CI95_tresh);

	upper95CI.push_back(rate_sample);
	sort(upper95CI.begin(), upper95CI.end());
	if(upper95CI.size()>CI95_tresh)
		upper95CI.erase(upper95CI.begin());

	lower99CI.push_back(rate_sample);
	sort(lower99CI.begin(), lower99CI.end());
	if(lower99CI.size()>CI99_tresh)
		lower99CI.resize(CI99_tresh);

	upper99CI.push_back(rate_sample);
	sort(upper99CI.begin(), upper99CI.end());
	if(upper99CI.size()>CI99_tresh)
		upper99CI.erase(upper99CI.begin());

}



