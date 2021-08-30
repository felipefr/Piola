#ifndef GPMATERIALS_H_
#define GPMATERIALS_H_

extern "C" {
	void __interfacefortran_MOD_executerelementc(int* , double*, double*, int* , 
						double*, int*, int*, int*, double*, double*, 
					 double*, double*, int*, double*, double*, double*);
}

extern "C" {
	int __determinant_MOD_lengthparam;
}

extern "C" {
	void __interfacefortran_MOD_executersymbolicc(int* , int*, double* ,
													int*, int*, int*);
}

inline void getLocalMatrix(int* id_Elem_Family, double* AE, double* BE, 
							int* MaxLRows, double* XLL, int* NDim, int* iDofT, 
							int* NodElt, double* Sol0, double* Sol1, 
					 double* CommonPar, double* Param, int* JParam, double* 
					 DelT, double* DTm, double* Time) {
	__interfacefortran_MOD_executerelementc(id_Elem_Family, AE, BE, MaxLRows, 
						XLL, NDim, iDofT, NodElt, Sol0, Sol1, 
						CommonPar, Param, JParam, DelT, DTm, Time);
}

inline void getLocalSymbolic(int* id_Elem_Family, int* Coupling, 
					double* CommonPar, int* iDofT, int* NDim, int* MaxLRows) {
	__interfacefortran_MOD_executersymbolicc(id_Elem_Family, Coupling,CommonPar,
													iDofT,NDim,MaxLRows);
}

#endif /* GPMATERIALS_H_ */
