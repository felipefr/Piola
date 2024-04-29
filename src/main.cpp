// main program of the   - Felipe Figueredo Rocha
#include <iostream>
#include "solver.hpp"
#include <armadillo>

//~ SET(ARMADILLO_DIR /home/felipe/PToolsV201311/armadillo-6.500.5/include)
//~ SET(ARMADILLO_LIB /home/felipe/PToolsV201311/armadillo-6.500.5/libarmadillo.so)

using namespace std;
using namespace arma;


int main(int argc, char **argv){	
	
	int n = 3;
	sp_mat A(n,n);
~ 
	A(0,0) = 1.0;
	A(1,1) = 2.0;
	A(2,1) = 4.0;
	A(2,2) = 3.0;
	 
	vec b(n);
	b(0) = 1.0;
	b(1) = 4.0;
	b(2) = 17.0;
	
	vec x = spsolve(A, b, "lapack");  // solve one system

	int n = 3;
	mat *A;
	A = new mat(n,n,fill::zeros);
~ 
	(*A)(0,0) = 1.0;
	(*A)(1,1) = 2.0;
	(*A)(2,1) = 4.0;
    (*A)(2,2) = 3.0;
	 
	vec *b;
	b = new vec(n,fill::zeros);
	 
	(*b)(0) = 1.0;
	(*b)(1) = 4.0;
	(*b)(2) = 17.0;
	 
	vec *x; 
	x = new vec(n,fill::zeros);
	
	//solve(*x, *A, *b);  // solve one system
	// cout << *x;


	Solver *s;
	
	s = new Solver;
	s->init();
	
	s->run();
	
	
}

