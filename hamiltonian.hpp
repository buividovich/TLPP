#ifndef _HAMILTONIAN_HPP_
#define _HAMILTONIAN_HPP_

using namespace std;

#include<random>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<fstream>
#include<algorithm>
#include<vector>
#include<limits>
#include<arpackpp/arrscomp.h>
#include <boost/program_options.hpp>

#include "linalg.hpp"
#include "square_lattice.hpp"
#include "ansi_io.hpp"
#include "timing.hpp"

namespace po = boost::program_options;

class Hamiltonian : public Square2DLattice
{
	private:
		int rng_engine_seed;
		std::ranlux48 rng_engine;
		std::uniform_real_distribution<double> rng_uniform_dist;
		std::normal_distribution<double> rng_normal_dist{0.0, 1.0};
	public:
		//Option holder
		Hamiltonian(int argc, char **argv);
		void init_parameters();
		void print_parameters();
		void get_suffix(char* cstr);
		po::options_description hamiltonian_options;
		string       CompoundName = "Generic";
		double       t[3];
		double       dt[3];
		double       lvec[3][2];
		
		double*      ts = NULL;  //Array of hoppings on all bonds
		/* Generating disorder */
		void         RandomizeHoppings();
		/* Hamiltonian definition */
		void         H(double* in, double* out);	 				//Hamiltonian
		void         JX(double* in, double* out);
		void         JY(double* in, double* out);
		void 	     Jmn2(double* res);
		void         Operator2Matrix(void (Hamiltonian::*O)(double* , double*), double* res);
		void         test_hermiticity(void (Hamiltonian::*O)(double* , double*), int sign=+1);
		/* Electric current definition */
		//void         J(double* in, double* out);	 				//Electric current, summed over all sites
		//Diagonalization routines
		void         diagonalise_H(bool noisy = false);
		void         test_eigensystem();
		double*      E    = NULL;									//Eigenenergies
		double*      psi  = NULL;									//The eigenvectors
};

#endif
