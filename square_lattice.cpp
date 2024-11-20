#include "square_lattice.hpp"

void Square2DLattice::idx2coords(uint idx, uint &x, uint &y)
{
	x = idx%L[0]; y = (idx/L[0])%L[1];
}

uint Square2DLattice::coords2idx(uint x, uint y)
{
	return x*B[0] + y*B[1];
}

Square2DLattice::Square2DLattice(int argc, char **argv)
{
	//Init options interface
	lattice_options.add_options()
		(          "L1",	   po::value<uint>(  &(L1)      )->default_value(    16), "L1")
		(          "L2",	   po::value<uint>(  &(L2)      )->default_value(    16), "L2");
		
	//Reading parameters from the command line
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(lattice_options).allow_unregistered().run(), vm);
	po::notify(vm);
	// Check if help is needed
	if(vm.count( "help" )){cout<<lattice_options<<endl;};
	
	init_parameters();
	print_parameters();
}

void Square2DLattice::init_parameters()
{
	//After reading the values from the command line, update all the constants	
	vol = L1*L2;
	
	L[0] = L1; L[1] = L2; L[2] =     L1; L[3] =     L2;
	B[0] =  1; B[1] = L1; B[2] =      1; B[3] =     L1;
	W[0] =  1; W[1] =  1; W[2] = (L1-1); W[3] = (L2-1);
}


void Square2DLattice::print_parameters()
{
	cout << ansi::magenta << "2D Lattice parameters: " << ansi::reset << endl;
	cout << "\t" << ansi::green << "Lattice size: \t" << ansi::yellow << L1  <<" x " << L2 << ansi::reset << endl;
}

void Square2DLattice::get_suffix(char* cstr)
{
	sprintf(cstr, "%ux%u", L1, L2);
}

uint Square2DLattice::shift(uint x, uint mu)
{
	uint tc   = (x/B[mu])%L[mu];
	uint nidx = x - tc*B[mu];
	tc    = (tc + W[mu])%L[mu];
	return (nidx + tc*B[mu]);
}

uint Square2DLattice::shift_by(uint x, uint mu, uint dx)
{
	uint tc   = (x/B[mu])%L[mu];
	uint nidx = x - tc*B[mu];
	tc    = (tc + dx)%L[mu];
	return (nidx + tc*B[mu]);
}

void Square2DLattice::test_shifts()
{
	for(uint x=0; x<vol; x++)
	{
		printf("x = (%02u, %02u)\t", (x%L[0]), (x/L[0])%L[1]);
		for(uint mu=0; mu<4; mu++)
		{
			uint xn = shift(x, mu);
			printf("[%01u->(%02u, %02u)] ", mu, (xn%L[0]), (xn/L[0])%L[1]);
		};
		printf("\n");
	};
}

uint Square2DLattice::link_index(uint x, uint mu)
{
	return ( ((mu==0)||(mu==1))? 2*x + mu : 2*shift(x, mu) + mu%2);
}

bool Square2DLattice::are_neighbours(uint x1, uint x2)
{
	if(x1==x2)
		return true;
	for(uint mu=0; mu<4; mu++)
		if(shift(x1,mu)==x2)
			return true;
	return false;
}


