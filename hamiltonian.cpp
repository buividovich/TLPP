#include "hamiltonian.hpp"

Hamiltonian::Hamiltonian(int argc, char **argv) : Square2DLattice(argc, argv)
{
	//Init options interface
	hamiltonian_options.add_options()
		(    "rng-seed",       po::value<int>(     &(rng_engine_seed)  )->default_value(       0), "Seed for the random number generator")
		(          "tA",	   po::value<double>(  &(t[0])             )->default_value(    57.0), "Hopping tA along direction 0 of the square lattice")
		(          "tB",	   po::value<double>(  &(t[1])      	   )->default_value(    57.0), "Hopping tB along direction 1 of the square lattice")
		(          "tC",	   po::value<double>(  &(t[2])             )->default_value(    57.0), "Diagonal hopping tC which turns square lattice to triangular one")
		(         "dtA",	   po::value<double>(  &(dt[0])            )->default_value(    28.0), "Dispersion of hopping tA")
		(         "dtB",	   po::value<double>(  &(dt[1])            )->default_value(    28.0), "Dispersion of hopping tB")
		(         "dtC",	   po::value<double>(  &(dt[2])            )->default_value(    28.0), "Dispersion of hopping tC")
		(        "eA_X",       po::value<double>(  &(lvec[0][0])       )->default_value(     0.1), "X component of the lattice vector that corresponds to A hopping")
		(        "eA_Y",       po::value<double>(  &(lvec[0][1])       )->default_value(     0.1), "Y component of the lattice vector that corresponds to A hopping")
		(        "eB_X",       po::value<double>(  &(lvec[1][0])       )->default_value(     0.1), "X component of the lattice vector that corresponds to B hopping")
		(        "eB_Y",       po::value<double>(  &(lvec[1][1])       )->default_value(    -0.1), "Y component of the lattice vector that corresponds to B hopping")
		("CompoundName",	   po::value<string>(  &CompoundName                                ), "Compound name (for output files)");
		
	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(hamiltonian_options).allow_unregistered().run(), vm);
	po::notify(vm);
	
	// Check if help is needed
	if(vm.count( "help" )){cout<<hamiltonian_options<<endl;};
	
	init_parameters();
	print_parameters();
}

void Hamiltonian::init_parameters()
{
	//Initialize the random number generator
	if(rng_engine_seed==0)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		int nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count();
		rng_engine.seed(nanos);
	}
	else
		rng_engine.seed(rng_engine_seed);
		
	//Initializing lattice vectors
	lvec[2][0] = lvec[0][0] - lvec[1][0];
	lvec[2][1] = lvec[0][1] - lvec[1][1];
	
	for(uint a=0; a<3; a++)
		dt[a] *= t[a];
	
	//Allocate memory for hoppings
	ts = new double[3*vol];
	RandomizeHoppings();
}

void Hamiltonian::print_parameters()
{
	cout << ansi::magenta << "Hamiltonian parameters for " << CompoundName << ansi::reset << endl;
	cout << "\t" << ansi::green << "Transfer integrals  tA,  tB,  tC (eV):                         \t" << ansi::yellow << t[0]  <<", " <<  t[1] << ", " <<  t[2] << ansi::reset << endl;
	cout << "\t" << ansi::green << "Their dispersions  dtA, dtB, dtC (eV):                         \t" << ansi::yellow << dt[0] <<", " << dt[1] << ", " << dt[2] << ansi::reset << endl;
	cout << "\t" << ansi::green << "Cartesian lattice vector for hopping tA (Angstrom):            \t" << ansi::yellow << "(" << lvec[0][0]  <<", " <<  lvec[0][1] << ")" << ansi::reset << endl;
	cout << "\t" << ansi::green << "Cartesian lattice vector for hopping tB (Angstrom):            \t" << ansi::yellow << "(" << lvec[1][0]  <<", " <<  lvec[1][1] << ")" << ansi::reset << endl;
	cout << "\t" << ansi::green << "Cartesian lattice vector for hopping tC (Angstrom):            \t" << ansi::yellow << "(" << lvec[2][0]  <<", " <<  lvec[2][1] << ")" << ansi::reset << endl;
}

void Hamiltonian::get_suffix(char* cstr)
{
	char lat_suffix[128];
	Square2DLattice::get_suffix(lat_suffix);
	sprintf(cstr, "%s_%s", CompoundName.c_str(), lat_suffix);
}

void Hamiltonian::RandomizeHoppings()
{
	#pragma omp parallel for
	for(uint x=0; x<vol; x++)
		for(uint a=0; a<3; a++)
			ts[3*x + a] = t[a] + dt[a]*rng_normal_dist(rng_engine);
}

void Hamiltonian::Operator2Matrix(void (Hamiltonian::*O)(double* , double*), double* res)
{
	double* B   = new double[vol];
	double* OB  = new double[vol];
	
	for(uint idx=0; idx<vol; idx++)
	{
		std::fill(B, B + vol, 0.0);
		B[idx] = 1.0;
		(this->*O)(B, OB);
		#pragma omp parallel for
		for(uint jdx=0; jdx<vol; jdx++)
			res[jdx*vol + idx] = OB[jdx];
	};
	
	delete [] OB;	delete [] B;
}

void Hamiltonian::H(double* in, double* out)
{
	std::fill(out, out + vol, 0.0);
	#pragma omp parallel for
	for(uint idx=0; idx<vol; idx++)
	{
		out[idx] -= ts[3*idx + 0]*in[shift(idx, 0)];
		out[idx] -= ts[3*idx + 1]*in[shift(idx, 1)];
		out[idx] -= ts[3*idx + 2]*in[shift(shift(idx, 0), 3)];
		
		out[idx] -= ts[3*shift(idx,2) + 0]*in[shift(idx, 2)];
		out[idx] -= ts[3*shift(idx,3) + 1]*in[shift(idx, 3)];
		
		uint jdx = shift(shift(idx, 1), 2);
		
		out[idx] -= ts[3*jdx + 2]*in[jdx];
	};
}

void Hamiltonian::test_hermiticity(void (Hamiltonian::*O)(double* , double*), int sign)
{
	double* M = new double[vol*vol];
	Operator2Matrix(O, M);
	double err = 0.0;
	for(uint idx=0; idx<vol; idx++)
		for(uint jdx=0; jdx<vol; jdx++)
			err = max(err, abs(M[idx*vol + jdx] - sign*M[jdx*vol + idx]));
	cout << ansi::magenta << "\tHermiticity error: " << ansi::yellow;
	printf("%2.2E", err);
	cout << ansi::reset << endl << flush;
	delete [] M;
}

void Hamiltonian::JX(double* in, double* out)
{
	std::fill(out, out + vol, 0.0);
	#pragma omp parallel for
	for(uint idx=0; idx<vol; idx++)
	{
		out[idx] -= lvec[0][0]*ts[3*idx + 0]*in[shift(idx, 0)];
		out[idx] -= lvec[1][0]*ts[3*idx + 1]*in[shift(idx, 1)];
		out[idx] -= lvec[2][0]*ts[3*idx + 2]*in[shift(shift(idx, 0), 3)];
		
		out[idx] += lvec[0][0]*ts[3*shift(idx,2) + 0]*in[shift(idx, 2)];
		out[idx] += lvec[1][0]*ts[3*shift(idx,3) + 1]*in[shift(idx, 3)];
		
		uint jdx = shift(shift(idx, 1), 2);
		
		out[idx] += lvec[2][0]*ts[3*jdx + 2]*in[jdx];
	};
}

void Hamiltonian::JY(double* in, double* out)
{
	std::fill(out, out + vol, 0.0);
	#pragma omp parallel for
	for(uint idx=0; idx<vol; idx++)
	{
		out[idx] -= lvec[0][1]*ts[3*idx + 0]*in[shift(idx, 0)];
		out[idx] -= lvec[1][1]*ts[3*idx + 1]*in[shift(idx, 1)];
		out[idx] -= lvec[2][1]*ts[3*idx + 2]*in[shift(shift(idx, 0), 3)];
		
		out[idx] += lvec[0][1]*ts[3*shift(idx,2) + 0]*in[shift(idx, 2)];
		out[idx] += lvec[1][1]*ts[3*shift(idx,3) + 1]*in[shift(idx, 3)];
		
		uint jdx = shift(shift(idx, 1), 2);
		
		out[idx] += lvec[2][1]*ts[3*jdx + 2]*in[jdx];
	};
}


void Hamiltonian::diagonalise_H(bool noisy)
{
	if(E==NULL    ){E     = new double[vol];};
	if(psi==NULL  ){psi   = new double[vol*vol];};

	TIMING_INIT;
	
	//Initialize storage for evecs
	Operator2Matrix(&(this->H), psi); //On input to LAPACKE_zheev, psi will contain the Hamiltonian matrix. Eigenvectors on output
	
	if(noisy) std::cout << "Running LAPACK_dsyev to find all " << vol << " lowest eigenstates " << endl << flush;
	
	TIMING_START;
	int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', vol, psi, vol, E);
	TIMING_FINISH;

	if(noisy) std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl; 
	if(res!=0)	cerr << ansi::red << "Something went wrong, LAPACKE_dsyev returned " << ansi::cyan << res << ansi::red << " !!!\n" << ansi::reset << endl << flush;

	if(noisy) std::cout << "Transposing the eigensystem... " << endl << flush;
	TIMING_START;
	for(uint i=0; i<vol; i++)
		for(uint j=i+1; j<vol; j++)
		{
			double tmp = psi[i*vol + j];
			psi[i*vol + j] = psi[j*vol + i];
			psi[j*vol + i] = tmp;
		};
	TIMING_FINISH;
	if(noisy) std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl; 
	
	if(noisy) std::cout << "Sorting the eigensystem... " << endl << flush;
	TIMING_START;
	sort_eigensystem(E, psi, vol, vol);
	TIMING_FINISH;
	if(noisy) std::cout << ansi::green << " ... Done in " << ansi::magenta << a_time << ansi::reset << " sec.\n" << endl;
}

void   Hamiltonian::test_eigensystem()
{
	double evec_err  = 0.0;
	double* tmp = new double[vol];
	for(uint iev=0; iev<vol; iev++)
	{
		H(psi + iev*vol, tmp);
		A_pluseq_bB(tmp, -E[iev], psi + iev*vol, vol);
		evec_err = max(evec_err, norm(tmp, vol));
	};
	delete [] tmp;
	
	cout << ansi::white << "Max. eigenvalue error is " << ansi::magenta << evec_err << ansi::reset << endl << flush;
	
	double ortho_err = orthonormality_norm(psi, vol, vol);

	cout << ansi::white << "Max. orthogonality error is " << ansi::magenta << ortho_err << ansi::reset << endl << flush;
}

void Hamiltonian::Jmn2(double* res)
{
	std::fill(res, res + vol*vol, 0.0);
	double* out  = new double[vol];
	for(uint n=0; n<vol; n++)
	{
		JX(psi + n*vol, out); // JX |n>
		for(uint m=0; m<vol; m++)
		{
			double Jmn = scalar_prod(psi + m*vol, out, vol); // <m| JX |n>
			res[m*vol + n] += Jmn*Jmn;
		};
		JY(psi + n*vol, out); // JY |n>
		for(uint m=0; m<vol; m++)
		{
			double Jmn = scalar_prod(psi + m*vol, out, vol); //<m| JY |n>
			res[m*vol + n] += Jmn*Jmn;
		};
	};
	delete [] out;
}