SRC = ./linalg.cpp ./hamiltonian.cpp ./square_lattice.cpp

HDR = $(SRC:.cpp=.hpp)
HDR += ./timing.hpp ./ansi_io.hpp 

CC = g++ -std=c++14 -O2 -fmax-errors=1 -fopenmp -I./ 
CC += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

LIB = -lm -lgfortran -lopenblas -lm
LIB +=  -lboost_program_options

CC += -I /opt/OpenBLAS/include/ -L /opt/OpenBLAS/lib/ 

tl: tl_main.cpp $(SRC) $(HDR)
	$(CC) $(SRC) ./$< $(LIB) -o ./tl
	
clean:
	rm -f -v ./tl