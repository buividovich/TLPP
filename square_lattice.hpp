#ifndef _SQUARE_LATTICE_HPP_
#define _SQUARE_LATTICE_HPP_

using namespace std;

#include <iostream>
#include <cstdio>
#include <boost/program_options.hpp>

#include "ansi_io.hpp"

namespace po = boost::program_options;

typedef unsigned int  uint; //Type for lattice coordinates and indexes

class Square2DLattice 
{ 
	public:
		Square2DLattice(){};
		//Option holder
		po::options_description lattice_options;
		//Actual parameters
		Square2DLattice(int argc, char **argv);
		void init_parameters();
		void print_parameters();
		void get_suffix(char* cstr);
    	uint L1 = 16;
    	uint L2 = 16;
    	uint vol;
    	uint W[4];
    	uint L[4];
		uint B[4];
    	uint shift(uint x, uint mu);
    	uint shift_by(uint x, uint mu, uint dx);
    	uint link_index(uint x, uint mu);
    	void test_shifts();
		bool are_neighbours(uint x1, uint x2);
		void idx2coords(uint idx, uint &x, uint &y);
		uint coords2idx(uint x, uint y);
};

#endif
