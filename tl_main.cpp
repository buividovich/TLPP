#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <vector>

#include "hamiltonian.hpp"
#include "timing.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
	TIMING_INIT;

	int      nthreads             = 0; //If 0, the number of threads will be set automatically
	string   datadir              = "./data/";
	//Physical temperature
	double   T                    = 25;
	uint    NT                    = 64;
	double   delta_relax          = 0.005;
	//Histogramming of the spectral function params
	double   wmax                 = 1.0; //Max. frequency for the histogram of the spectral functions
	double   dw                   = 0.01;  //Frequency resolution for the histogram
	char     spf_hist_fname[512]  = "";
	char           GE_fname[512]  = "";
	//Suffixes for file naming
	char	 suffix[512]          = "";
	uint     nconfigs             = 1;
	
	po::options_description general_options("Algorithm options");
	general_options.add_options()
		("help,h", "produce this help message")
		("nthreads",	   po::value<int>(        &(nthreads))->default_value(        0), "Number of OpenMP threads to use, 0 = automatic choice"        )
		("datadir", 	   po::value<string>(                   &datadir               ), "Directory for data output"                                    )
		("T",              po::value<double>(            &(T))->default_value(    0.025), "Temperature"                                                  )
		("NT",             po::value<uint>(             &(NT))->default_value(       64), "Number of discrete values of imaginary time between 0 and NT" )
		("nconfigs",       po::value<uint>  (     &(nconfigs))->default_value(        1), "Number of disorder configurations to average over"            )
		("wmax", 	       po::value<double>(         &(wmax))->default_value(      1.0), "Max. frequency for spectral function histogramming"           )
		("dw", 	           po::value<double>(           &(dw))->default_value(    0.001), "Bin width for spectral function histogramming"                )
		("delta-relax",    po::value<double>(  &(delta_relax))->default_value(    0.005), "Relaxation time for RTA mobility calculation"                 );
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(general_options).allow_unregistered().run(), vm);
	po::notify(vm);
	
	// Check if help is needed
	if(vm.count( "help" )){cout<<general_options<<endl; return 1;};
	
	//Printing general options
	double beta = 1.0/T;
	
	nthreads = (nthreads==0? omp_get_max_threads() : nthreads);
	omp_set_num_threads(nthreads);
	cout << endl << ansi::cyan << "\t This code implements SSH Hamiltonian with static bond disorder " << endl << endl << flush;
	cout << ansi::white << "Using " << ansi::magenta << nthreads << ansi::white << " OpenMP threads" << ansi::reset << endl;
	
	cout << ansi::green << "Directory for data output:                                 \t" << ansi::white   << datadir               << ansi::reset << endl;
	cout << ansi::green << "Temperature (eV):                                          \t" << ansi::magenta << T                                           ;
	cout << ansi::cyan  << " (beta = " << ansi::magenta << beta << ansi::cyan << ")"                                          << ansi::reset << endl; 
	cout << ansi::green << "Inverse relaxation time for RTA mobility calculation (eV): \t" << ansi::magenta << delta_relax           << ansi::reset << endl;
	cout << ansi::green << "No. of Euclidean time steps:                               \t" << ansi::magenta << NT                    << ansi::reset << endl;
	cout << ansi::green << "Max. frequency for spectral function histogramming (eV):   \t" << ansi::magenta << wmax                  << ansi::reset << endl;
	cout << ansi::green << "Bin width for spectral function histogramming (eV):        \t" << ansi::magenta << dw                    << ansi::reset << endl;
	
	//Diagonalizing the spin chain without any magnetic field ...
	Hamiltonian* OSH = new Hamiltonian(argc, argv);
	
	char h_suffix[256]; OSH->get_suffix(h_suffix);
	sprintf(suffix, "%s_T%2.4lf", h_suffix, T);
	
	sprintf(spf_hist_fname, "%s/spf_hist_static_%s.dat", datadir.c_str(), suffix);
	sprintf(      GE_fname, "%s/GE_static_%s.dat",       datadir.c_str(), suffix);
	
	cout << ansi::green << "Suffix for output files:                  \t" << ansi::white   << suffix                << ansi::reset << endl;
	cout << ansi::green << "Output file for the SPF histogram:        \t" << ansi::white   << spf_hist_fname        << ansi::reset << endl;
	cout << ansi::green << "Output file for GE(tau):                  \t" << ansi::white   << GE_fname              << ansi::reset << endl;
	cout << endl;
	
	OSH->test_hermiticity(&(OSH->H));
	OSH->test_hermiticity(&(OSH->JX), -1);
	OSH->test_hermiticity(&(OSH->JY), -1);
	
	double* Jmn2 = new double[OSH->vol*OSH->vol];
	double* GE   = new double[NT+1];
	std::fill(GE, GE + NT + 1, 0.0);
	
	uint NW = (uint)ceil(wmax/dw) + 1;
	double* SPF  = new double[NW];
	std::fill(SPF, SPF + NW, 0.0);
	
	double  Z    = 0.0;  

	double  mobility_rta = 0.0;
	
	FILE* GE_file = fopen(GE_fname, "w");
	if(GE_file==NULL) cout << ansi::red << "Could not open the file " << ansi::cyan << GE_fname << ansi::red << " for writing" << ansi::reset << endl << flush;
	
	for(uint iconf=0; iconf<nconfigs; iconf++)
	{
		OSH->RandomizeHoppings();
		OSH->diagonalise_H();
		//OSH->test_eigensystem();
		OSH->Jmn2(Jmn2);
		
		for(uint m=0; m<OSH->vol; m++)
			for(uint n=0; n<OSH->vol; n++)
			{
				double w = OSH->E[n] - OSH->E[m];
				
				if(w>=0.0)
				{
					uint ibin = (uint)floor(w/dw);
					if(ibin<NW) SPF[ibin] += Jmn2[m*OSH->vol + n]*(exp(-beta*OSH->E[m]) + exp(-beta*OSH->E[n]));
				};
				
				for(uint it=0; it<=NT; it++)
				{
					double tau = (double)it*beta/(double)NT;
					GE[it] += Jmn2[m*OSH->vol + n]*exp(-tau*OSH->E[m] - (beta - tau)*OSH->E[n]);
				};

				mobility_rta += Jmn2[m*OSH->vol + n]*exp(-beta*OSH->E[n])/(delta_relax*delta_relax + w*w);
			};
				
						
		for(uint n=0; n<OSH->vol; n++)
			Z += exp(-beta*OSH->E[n]);
			
		cout << "." << flush;
	};
	cout << endl;

	//This implements Eq. (1) in DOI:10.1038/NMAT4970
	mobility_rta = mobility_rta*delta_relax/(2.0*Z*T);
	
	cout << ansi::cyan << "Partition function Z: " << ansi::yellow << Z << ansi::reset << endl;
	cout << ansi::cyan << "RTA mobility:         " << ansi::yellow << mobility_rta/6.58 << " cm^2/(V*s)" << ansi::reset << endl;
	
	cout << endl << ansi::cyan << "GE(tau): " << ansi::reset << endl;
	for(uint it=0; it<=NT; it++)
	{
		double tau = (double)it*beta/(double)NT;
		GE[it] /= (Z*OSH->vol);
		if(GE_file!=NULL) fprintf(GE_file, "%2.4lf\t%+2.6E\n", tau, GE[it]);
		printf("%2.4lf\t%+2.6E\n", tau, GE[it]);
	};
	cout << endl << endl;
	
	if(GE_file!=NULL) fclose(GE_file);
	
	FILE* SPF_file = fopen(spf_hist_fname, "w");
	if(SPF_file==NULL) cout << ansi::red << "Could not open the file " << ansi::cyan << spf_hist_fname << ansi::red << " for writing" << ansi::reset << endl << flush;
	if(SPF_file!=NULL)
	{
		for(uint iw=0; iw<NW; iw++)
		{
			double w = ((double)iw + 0.5)*dw;
			SPF[iw] /= (Z*((double)OSH->vol)*dw);
			fprintf(SPF_file, "%2.4lf\t%+2.6E\n", w, SPF[iw]);
		};
		fclose(SPF_file);
	};
	
	return EXIT_SUCCESS;
}