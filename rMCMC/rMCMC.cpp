/*
 * rMCMC.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: claudebherer
 */

#include <ctime>
#include <iostream>
#include <limits>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>  // shuffle

#include "readFile.h"
#include "interval.h"
#include "event.h"
#include "set_intervals.h"
#include "LogFactorial.h"
#include "read_nbmeioses.h"

using namespace std;

double interval::prior_recomb_rate = 1;
double interval::prior_variance_rate = 1;
int interval::CI95_tresh = 25;
int interval::CI99_tresh = 10;


std::default_random_engine interval::generator(1);
std::default_random_engine event::generator(2);

void computeMLE(vector<event> &events, vector<interval> &intervals, double &mean_MLE)
{
	// Compute the mean MLE of recombination rate across whole region
	// 1. get proportional contribution of each event to each interval
	vector<double> fraction(intervals.size(),0);
	for (int i=0; i<events.size(); i++)
	{
		double total_d = events[i].max_interval_pos-events[i].min_interval_pos;

		for (unsigned int ui=0; ui<events[i].N_intervals; ui++)
		{
			double d_bp = (intervals[events[i].min_interval_idx+ui].stoppos - intervals[events[i].min_interval_idx+ui].startpos)/1.0;
			fraction[events[i].min_interval_idx+ui] += d_bp/total_d;

		}
	}
	// 2. Calculate the MLE of recombination rate for each interval
	mean_MLE = 0.0;
	for (int i=0; i<intervals.size(); i++)
	{
		double rfr = fraction[i]/intervals[i].M;
		double MLE = 25.0 * log((1.0+2.0*rfr) / (1.0-2.0*rfr));	// cM
		mean_MLE += MLE;
		MLE /= ((intervals[i].stoppos - intervals[i].startpos) / 1000000.0);	// cM/Mb
		intervals[i].set_MLE(MLE);

		//cerr << "interval " << i << " MLE " << intervals[i].MLE << endl;
	}
	mean_MLE /= ((intervals[intervals.size()-1].stoppos - intervals[0].startpos) / 1000000.0);
}


int main(int argc, char* argv[])
{
	// time check
	time_t start_time, end_time;
	time(&start_time);
	struct tm * timeinfo;
	timeinfo = localtime(&start_time);
	cerr << endl << "Start time: " << asctime(timeinfo) << endl;

    string filename=""; // filename (and path) to the events file
    string out_prefix = "out";
    int M=-1;			// nb of meioses
    int nchains = 100000;	// nb of iterations of the MCMC chain
    int burnin = 33000;     // nb of burn-in iterations
	double prior_rate = 1.0; // prior recombination rate of 1 cM/Mb
    double prior_var = 5.0; // prior variance set by default to 5 (var=0.66 at 1Mb scale in HapMap2)
    int samp = 1000;
    int seed = (int)start_time;
	double HDI_region = 0.95;
	bool prior_rate_set = false;
	string M_filename=""; // Path/filename to the nb of meioses file

    if (argc < 2)
    {
    	cerr << "Arguments of rMCMC:" << endl;
    	cerr <<	"-i /path/filemame (event file)" << endl;
    	cerr << "-m nb of meioses" << endl << "-c nb of iterations of the MCMC chain" << endl;
    	cerr << "-r prior mean (cM/Mb)" << endl;
    	cerr << "-v prior variance at 1Mb scale (cM^2)" << endl;
    	cerr << "-c number of iterations" << endl;
    	cerr << "-b nb of burn-in iterations" << endl;
    	cerr << "-d iterations between samples" << endl;
    	cerr << "-s random seed" << endl;
    	cerr << "-h HDI region probability" << endl;
    	cerr << "-o output_prefix" << endl;
    	cerr << "-nbmeioses /path/filemame (nb meioses file)" << endl;
    }

	unsigned int ui=1;
    while (ui+1 < argc)
    {
    	string in_str = argv[ui];
   	    string in_str2;
    	if (in_str == "-i") // argument 1 : input file : event file
    	{	// Input file
    		in_str2 = argv[ui+1];
    		filename = in_str2.c_str();
    	}
    	else if (in_str == "-m") // argument 2 : total nb of meioses
    	{	// Number of meioses
    		in_str2 = argv[ui+1];
    		M = atoi(in_str2.c_str());
    	}
    	else if (in_str == "-c") // argument 3 : nb of chains (iterations) of the gibbs samplers
        {  	// Number of iterations
        	in_str2 = argv[ui+1];
        	nchains = atoi(in_str2.c_str());
        }
    	else if (in_str == "-b") // argument 4 : nb of burn-in iterations
        {  	// Burn in 
        	in_str2 = argv[ui+1];
        	burnin = atoi(in_str2.c_str());
        }
     	else if (in_str == "-v") // argument 5 : prior variance
    	{	// Variance
    	   	in_str2 = argv[ui+1];
    	   	prior_var = atof(in_str2.c_str());
    	}
    	else if (in_str == "-r")
    	{	// Prior mean rate
    		in_str2 = argv[ui+1];
        	prior_rate = atof(in_str2.c_str()); // argument 6 : prior mean rate
        	prior_rate_set = true;
    	}
    	else if (in_str == "-s")
    	{	// Read in seed.
    		in_str2 = argv[ui+1];
        	seed = atoi(in_str2.c_str());
    	}
    	else if (in_str == "-d")
    	{	// Read in sample rate.
    		in_str2 = argv[ui+1];
        	samp = atoi(in_str2.c_str());
    	}
    	else if (in_str == "-o")
    	{
	    	out_prefix = argv[ui+1];
    	}
    	else if (in_str == "-h")
    	{	// Read in HDI region probability
    		in_str2 = argv[ui+1];
        	HDI_region = atof(in_str2.c_str());
    	}
    	else if (in_str == "-nbmeioses")
    	{	// Read in /path/filemame (nb meioses file)
       		in_str2 = argv[ui+1];
        	M_filename = in_str2.c_str();
    	}
    	ui++;
    }
    
    if ((filename == "") || (M == -1))
    {
    	cerr << "Error: missing parameters." << endl; exit(1);
	}

    // Read events file and construct event objects
    vector<event> events;
  	readFile(filename, events);

  	unsigned int nb_events = (unsigned int)events.size();
  	cerr << "Number of events: " << nb_events << endl;

 	// Read nbmeioses file and build intervals
 	vector<interval> intervals;
  	read_nbmeioses(M_filename, intervals);
	cerr << "nbmeioses file read" << endl;
 	set_intervals(events, intervals); // find overlap btw intervals and events

 	unsigned int nb_intervals = (unsigned int)intervals.size();
	cerr << "Number of sequence intervals: " << nb_intervals << endl;

	/*
	// Print to check
	cerr << "int i\tstart\tstop\tnb meioses\n";
	for(int i=0; i<20; i++)
	{
		cerr << i << "\t" << intervals[i].startpos << "\t" << intervals[i].stoppos << "\t" << intervals[i].M << endl;
	}
	*/

	event::set_seed(seed);
 	interval::prior_recomb_rate = prior_rate;
 	interval::prior_variance_rate = prior_var;
 	interval::set_seed(seed+1);

	cerr << "Filename: " << filename << endl;
	cerr << "Output Prefix: " << out_prefix << endl;
    cerr << "Prior mean: " << prior_rate << " cM/Mb" << endl;
    cerr << "Prior variance at 1Mb scale: " << prior_var << " cM^2" << endl;
    cerr << "Seed: " << seed << endl;
    cerr << "Meioses: " << M << endl;
    cerr << "Iterations: " << nchains << endl;
    cerr << "Burn: " << burnin << endl;
    cerr << "Sample: " << samp << endl;
    cerr << "HDI region probability: " << HDI_region << endl;

	if (prior_rate_set == false) // if prior_rate is not passed in command line
	{
		double mean_MLE;
		computeMLE(events, intervals, mean_MLE);
		cerr << "Mean MLE: " << mean_MLE << " cM/Mb" << endl;
		prior_rate = mean_MLE;
		interval::prior_recomb_rate = prior_rate;
	}


 	// Starting state for the MCMC: Set Gamma priors for each sequence interval and draw lambda
 	for (int i=0; i<nb_intervals; i++)
 	{
		intervals[i].set_gamma_priors();
		intervals[i].lambdaij = intervals[i].lambda0;
		intervals[i].set_alpha_and_beta(intervals[i].alpha0, intervals[i].beta0);
	}
	
	cerr << "Region range: " << intervals[0].startpos << "-" << intervals[intervals.size()-1].stoppos << " bp" << endl;

 	// Gibbs sampler (nchains)
	// Step 1. Sample the location of recombination events given lambdas
	// Step 2. Sample the rec. rates in each interval

 	// Output vectors
  	vector<int> sample_count(intervals.size(),0);
	vector<double> sample_sum(intervals.size(),0);
  	vector<double> sample_sumsq(intervals.size(),0);

  	// Non-parametric 95%CI point estimates
  	interval::CI95_tresh = ((nchains-burnin)/samp) * 0.025;
  	cout << "Non-parametric 95%CI point estimates: " << interval::CI95_tresh << endl;
  	interval::CI99_tresh = ((nchains-burnin)/samp) * 0.005;
  	cout << "Non-parametric 99%CI point estimates: " << interval::CI99_tresh << endl;

	cerr << endl << "Iteration\tLog_LK\tMean_rate" << endl;
  	for (int nchain=0; nchain<nchains; nchain++)
  	{
  		// Step 1. Sample the location of recombination events given lambdas P(rec|lambdai/sum(lambdas))
  		// 1.1. Draw events and update the vector of the nb of events currently in interval i
  		vector<int> vec_nb_events(intervals.size(), 0); // vector of the nb of events currently in interval i
  	 	for (int i=0; i<nb_events; i++)
  	 	{
  	 		events[i].draw_new_current_interval(intervals);
  	 		vec_nb_events[events[i].current_interval_idx]++;
  	 	}

 	 	// 1.2. Update alphaij, betaij
  	 	for (int i=0; i<nb_intervals; i++)
  	 	{
  	 		intervals[i].alphaij = intervals[i].alpha0 + vec_nb_events[i];
  	 		// Not sure if 2 is needed here. 
  	 		intervals[i].betaij = intervals[i].beta0 + intervals[i].M;

	  	 	// Step 2. Sample the rec. rates in each interval
   			intervals[i].draw_new_lamba();
  		}

	  	// Take samples every 1000 iterations, after a given nb of burn-in samples
	  	if (((nchain+1) > burnin) && (((nchain+1) % samp) == 0))
	  	{
			// Compute the log likelihood of rates to test for convergence
			double log_likelihood = 0.0;
		  	for (int i=0; i<nb_intervals; i++)
		  	{
		  		// Calculate probability mass function
		  		if (intervals[i].lambdaij > 0)
		  		{
		  			double param = intervals[i].M*intervals[i].lambdaij;	// Mean of Poisson distribution
		  			log_likelihood += ((vec_nb_events[i]*log(param)) - param) - LogFactorial(vec_nb_events[i]); // probability mass function
		  		}
		  		else 
		  			log_likelihood += -9999999.9;
			}

		  	double mean_rate = 0.0;		// Mean rate across whole region
	  		for (int i=0; i<nb_intervals; i++)
	  		{
	  			sample_count[i]++;
	  			double rate = 25.0 * log((1.0+2.0*intervals[i].lambdaij) / (1.0-2.0*intervals[i].lambdaij));	// cM
	  			mean_rate += rate;
	  			rate /= ((intervals[i].stoppos - intervals[i].startpos) / 1000000.0);	// cM/Mb
	  			sample_sum[i] += rate;
	  			sample_sumsq[i] += (rate * rate);


	  			intervals[i].set_CI(rate);

			}
			mean_rate /= ((intervals[intervals.size()-1].stoppos - intervals[0].startpos) / 1000000.0);
	  		cerr << (nchain+1) << "\t" << log_likelihood << "\t" << mean_rate << endl;
	  		
	  		for (int i=0; i<nb_events; i++)
	  			events[i].take_interval_sample();

	  	} // fin if sample

  	} // fin boucle nchain (Gibbs sampler)

	// Output Map estimate
	string out_filename1 = out_prefix + "-rates.txt";
	ofstream out(out_filename1);
	// Output a table with: interval_ID, interval_startpos, interval_stoppos, rfr,
  	// nb of samples (count of lambdas values sampled), sum of lambda values and sum of squares of lambda values
  	out << "#int_ID\tstart_pos\tstop_pos\tprior_recomb_fraction\tmean\tvar\tlower_95\tupper_95\tlower_99\tupper_99" << endl;
  	for (int i=0; i<sample_count.size();i++)
  	{
  		out << i << "\t";
  		out << intervals[i].startpos << "\t";
  		out << intervals[i].stoppos << "\t";
  		out << intervals[i].r0 << "\t";
  		double mean = sample_sum[i] / sample_count[i];
  		double var = (sample_count[i]*sample_sumsq[i] - (sample_sum[i]*sample_sum[i])) / (sample_count[i] * (sample_count[i]-1));
  		//double sd = sqrt(var);
  		double l95 = intervals[i].lower95CI[interval::CI95_tresh-1];
  		double u95 = intervals[i].upper95CI[0];
  		double l99 = intervals[i].lower99CI[interval::CI99_tresh-1];
  		double u99 = intervals[i].upper99CI[0];
  		out << mean << "\t" << var << "\t" << l95 << "\t" << u95 << "\t" << l99 << "\t" << u99 << endl;
  	}
  	out.close();
  	
  	// Output Event posterior estimates

  	string out_filename2 = out_prefix + "-events.txt";
  	ofstream out2(out_filename2);
  	out2 << "#Event_ID\tOrig_Start\tOrig_End\tOrig_Size\tLow_95%_HDI\tHigh_95%_HDI\tHDI_Size" << endl;
  	for (int i=0; i<nb_events; i++)
  	{
  		double sum = accumulate(events[i].interval_samples.begin(), events[i].interval_samples.end(), 0);
  		vector<double> post_cdf(events[i].interval_samples.size()+1, 0.0);
  		unsigned int LHDI = events[i].min_interval_pos;
  		unsigned int UHDI = events[i].max_interval_pos;
  		double thresh1 = (1.0 - HDI_region) * 0.5;
  		double thresh2 = 1.0 - thresh1;
  		for (unsigned int ui=0; ui<events[i].interval_samples.size(); ui++)
  		{
  			post_cdf[ui+1] = post_cdf[ui] + events[i].interval_samples[ui] / sum;
  			if ((post_cdf[ui] <= thresh1) && (post_cdf[ui+1] >= thresh1))
  			{	// Lower limit
  				unsigned int pos_idx = events[i].min_interval_idx + ui;
  				double pos1 = intervals[pos_idx].startpos;
  				double pos2 = intervals[pos_idx].stoppos;  				
  				LHDI = (unsigned int) ((pos1 + ((thresh1 - post_cdf[ui]) / (post_cdf[ui+1] - post_cdf[ui])) * (pos2 - pos1)) + 0.5);
  			}
  			else if ((post_cdf[ui] <= thresh2) && (post_cdf[ui+1] >= thresh2))
  			{	// Upper limit
  				unsigned int pos_idx = events[i].min_interval_idx + ui;
  				double pos1 = intervals[pos_idx].startpos;
  				double pos2 = intervals[pos_idx].stoppos;  				
  				UHDI = (unsigned int) ((pos1 + ((thresh2 - post_cdf[ui]) / (post_cdf[ui+1] - post_cdf[ui])) * (pos2 - pos1)) + 0.5);  			
  			}
  		}
  		
  		out2 << i << "\t" << events[i].min_interval_pos << "\t" << events[i].max_interval_pos;
  		out2 << "\t" << events[i].max_interval_pos - events[i].min_interval_pos << "\t";
  		out2 << LHDI << "\t" << UHDI << "\t" << UHDI - LHDI << endl;
  	}
  	out2.close();

   	time(&end_time);
   	time_t seconds;
   	seconds = difftime(end_time, start_time);

   	// Convert time
   	unsigned int minutes, hours, secs_left, mins_left;
   	hours = (unsigned int)seconds / 3600;
   	minutes = (unsigned int)seconds / 60;
   	mins_left = minutes % 60;
   	secs_left = seconds % 60;

	timeinfo = localtime(&end_time);
	cerr << endl << "End time: " << asctime(timeinfo) << endl;
   	cerr << "Total run time: " << hours << ":" << mins_left << ":" << secs_left << endl;

	return 0;
}
