#ifndef _COHERENT_STATES_HPP_
#define _COHERENT_STATES_HPP_

using namespace std;

#include "hamiltonian.hpp"

class CSH : public Hamiltonian
{
	public:
		//Option holder
		CSH(int argc, char **argv);
		void init_parameters();
		void print_parameters();
		void get_suffix(char* cstr);
		po::options_description csh_options;
		/* Generating disorder */
		//void         RandomizeHoppings();
		/* Hamiltonian definition */
		//void         H(double* in, double* out);	 				//Hamiltonian
		//void         JX(double* in, double* out);
		//void         JY(double* in, double* out);
		//void 	     Jmn2(double* res);
		//void         Operator2Matrix(void (Hamiltonian::*O)(double* , double*), double* res);
		//void         test_hermiticity(void (Hamiltonian::*O)(double* , double*), int sign=+1);
		/* Electric current definition */
		//void         J(double* in, double* out);	 				//Electric current, summed over all sites
		//Diagonalization routines
		//void         diagonalise_H(bool noisy = false);
		//void         test_eigensystem();
		//double*      E    = NULL;									//Eigenenergies
		//double*      psi  = NULL;									//The eigenvectors
};

#endif
