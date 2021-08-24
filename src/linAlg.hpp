// Implements non-natives linear algebra classes (Mat and Vec) - Felipe Figueredo Rocha
#ifndef _linalg_hpp
#define _linalg_hpp

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>  
#include <sstream>
#include <cmath>
//~ #include <gsl/gsl_matrix.h>
//~ #include <gsl/gsl_vector.h>
//~ #include <gsl/gsl_blas.h>
//~ #include <gsl/gsl_linalg.h>
#include <armadillo>

//~ using namespace arma;
using namespace std;

#define PRECISAO 12 


// defines shortcuts for some most useds object types
template <class T> class Vec;
template <class T> class Mat;
typedef double SGPreal; // long double is not supported by GSL
typedef Vec<SGPreal> SGPrealVec;
typedef Mat<SGPreal> SGPrealMat;
typedef Vec<int> intVec;
typedef Mat<int> intMat;

// class numerical vector
template <class T> class Vec{
	public:
	T *v;
	int n;
	Vec() { n = 0; }
	Vec(int nn) { n= nn; v = new T[n]; (*this)=0.0; }	
	void readNumberBlock(ifstream &file); // reads a block of n (size of the vector) numbers in a file
	void writeNumberBlock(ofstream &file); // idem to the last, but writes
	void writeNumberBlockRect(ofstream &file, int m);
	void print();
	void printH();
	T max();
	T min();
	T sum();
	T amax();
	T& operator()(int i) {return v[i];}
	T& operator[](int i) {return v[i];}
	void operator=(T a) { for(int i=0; i<n ; i++) v[i]=a;}
	void operator=(Vec &w) { for(int i=0; i<n ; i++) v[i]=w[i];}
	void operator=(arma::vec &w) { for(int i=0; i<n ; i++) v[i]=w[i];}
};

template<class T> void Vec<T>:: print(){	
	for(int i = 0; i<n ; i++) cout << v[i] << endl;
}

template<class T> void Vec<T>:: printH(){	
	for(int i = 0; i<n ; i++) cout << v[i] << " ";
	cout << endl;
}


template <class T> void Vec<T> :: readNumberBlock(ifstream &file){	
	//~ file << fixed << setprecision(PRECISAO);
	for(int i = 0; i<n ; i++) file  >> v[i] ;	
}

template <class T> void Vec<T> :: writeNumberBlock(ofstream &file){	
	file << scientific << setprecision(PRECISAO);
	for(int i = 0; i<n ; i++) file  << v[i] << endl ;	
}

template <class T> void Vec<T> :: writeNumberBlockRect(ofstream &file,int m){	
	if(n%m==0){
		int nn = n/m;
		int ip = 0;
		file << scientific << setprecision(PRECISAO);
		for(int i = 0; i<nn ; i++){ 
			ip = i*m;
			for(int j = 0; j<m ; j++){
				file  << v[ip + j] << " " ;
			}
			file << endl ;
		}
	} else {
		writeNumberBlock(file);
	}
}	

template <class T> T Vec<T> :: max(){	
	T vmax = -9999.0;
	
	for(int i = 0; i<n ; i++){
		if(v[i]>vmax) vmax = v[i];
	}	
	
	return vmax;
}

template <class T> T Vec<T> :: sum(){	
	T vsum = 0.0;
	
	for(int i = 0; i<n ; i++) vsum += v[i];
	
	return vsum;
}

template <class T> T Vec<T> :: amax(){	
	T vmax = 0.0;
	T aux;
	
	for(int i = 0; i<n ; i++){
		aux = abs(v[i]);
		if(aux>vmax) vmax = aux;
	}	
	
	return vmax;
}


// class numerical matrix as a spelization of the Vec. All the elements are in a Vec, but can be acessed with two indices. Row-major convention (C convention)
template <class T> class Mat : public Vec<T>{
	public:
	int m1,m2; // m1 is Nrows, m2 is Ncolumns 
	Mat() {}
	Mat(int mm1, int mm2):Vec<T>(mm1*mm2) { m1 = mm1; m2 = mm2;  }
	Mat(Vec<T> &V, int mm1, int mm2) { m1 = mm1; m2 = mm2; this->v = V.v; }
	
	//~ T& operator()(int i,int j) {return this->v[j*m1 + i];} Fortran style
	T& operator()(int i,int j) {return this->v[i*m2 + j];} // C style, compatibility with GSL
	void operator=(T a) { for(int i=0; i<this->n ; i++) this->v[i]=a;}
	void operator=(Mat &w) { for(int i=0; i<this->n ; i++) this->v[i]=w->v[i];}
	void prettyPrint();
	void prettyPrintClean();
	void solve(SGPrealVec &x,SGPrealVec &b); // solves linear system using GSL
};

template<class T> void Mat<T> :: prettyPrint(){
	//~ cout << m1 << " " << m2 << endl;
	for (int i=0;i<m1;i++){ 
        for (int j=0;j<m2;j++) cout << "a[" << i << "][" << j << "]=" << (*this)(i,j) << "  " ;
        cout << endl; 
    }
}

template<class T> void Mat<T> :: prettyPrintClean(){
	//~ cout << m1 << " " << m2 << endl;
	for (int i=0;i<m1;i++){ 
        for (int j=0;j<m2;j++) cout << (*this)(i,j) << "  " ;
        cout << endl; 
    }
}

template<class T> T dot_product(Vec<T> &u,Vec<T> &w){
	T dot = 0.0;
	
	for(int i=0; i<u.n ; i++) dot += u(i)*w(i);
	
	return dot;
}

template<class T> void Mat<T> :: solve(SGPrealVec &x,SGPrealVec &b){ 
	// creates GSL objects
	//~ gsl_matrix_view gslA=gsl_matrix_view_array(this->v,this->m1,this->m2); 
	//~ gsl_vector_view gslB=gsl_vector_view_array(b.v,b.n);
	//~ gsl_vector_view gslX=gsl_vector_view_array(x.v,x.n);
	//~ gsl_permutation *gslP=gsl_permutation_alloc(b.n);
	//~ 
	//~ int sig;
      //~ 
    //~ // solve system
	//~ gsl_linalg_LU_decomp(&gslA.matrix,gslP,&sig);
    //~ gsl_linalg_LU_solve(&gslA.matrix,gslP,&gslB.vector,&gslX.vector);
}

double computeError(arma::vec &v, SGPrealVec &v0){
	double error = 0.0;	
	arma::vec e(v0.n);
	
	for(int i=0;i<v0.n;i++) e(i) = v(i) - v0(i);
	
	double emax = arma::abs(e).max();
	double v0max = v0.amax();
	
	if(v0max >0.0){
		error = emax/v0max;
	}
	else{
		error = emax;
	}
	
	return error;
}

double computeError(SGPrealVec &v,SGPrealVec &v0){
	
	SGPrealVec e(v.n);
	double error = 0.0;	
	
	for(int i=0; i<v.n ; i++) e(i) = v(i) - v0(i);
	
	double emax = e.amax();
	double v0max = v0.amax();
	
	if(v0max >0.0){
		error = emax/v0max;
	}
	else{
		error = emax;
	}
	
	return error;
}


#endif
