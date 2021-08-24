// implements the solver class - Felipe Figueredo Rocha 
#ifndef _solver_hpp
#define _solver_hpp

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>  
#include <sstream>
#include <cmath>
#include "utils.hpp"
#include "linAlg.hpp"
#include "data.hpp"
#include "simulSettings.hpp"
#include "solvergp.h"
#include <armadillo>

//~ #define DEBUG if(false)
#define DEBUG if(true)

using namespace arma;
using namespace std;

class Solver{
	public:
	// see constructor for further explanation about the variabls
	Data *d;
	SimulSettings *s;
	mat *A; // stiffness matrix
	vec *b; // force vector
	vec *Sol; // force vector
	vector<intVec*> *Coupling; // Coupling vector, works similar to flagDirich
	Solver() {}
	void init();
	void prepareMatrices(); // alloc linear system matrices and vector
	void createCoupling(); 
	void linsolve(); // call the routine of solving the lin System
	void enforceDirichlet(SGPrealVec *vDirich, intVec *flagDirich, double pen); // modify matrices to enforce Dirichlet
	void assembly(SGPrealMat *A_L, SGPrealVec *b_L, int e); // assembly local in global matrices
	void assemblyCoupling(intVec *Cvec, intMat *Cmat_L, int e);
	void globalToLocal(SGPrealVec *XLL,SGPrealVec *Sol0E, SGPrealVec *Sol1E,SGPrealVec *ParamE, int e); // fills the e^{th} of the mesh in a auxiliary object
	void localToGlobal(SGPrealVec *ParamE,int e); // fills the e^{th} of the mesh in a auxiliary object
	void run();
};

void Solver:: init(){
	s = new SimulSettings;
	s->loadSettings();
	d = new Data;
	d->loadSGPmesh();
	d->loadSGPdirich(s->NSubsteps);
	d->loadSGPinifile();
	d->loadSGPparam();
	prepareMatrices();
	createCoupling();
}

void Solver:: prepareMatrices(){
	int n = d->Ndof*d->Nnodes;

 	A = new mat(n,n);
	b = new vec(n);
	Sol = new vec(n);

}

void Solver:: createCoupling(){
	int n = d->Ndof*d->Nnodes;
	int NodElt = d->eNNod->max();
	int MaxLRows = NodElt*d->Ndof;
	intMat *CouplingE = new intMat(MaxLRows,MaxLRows);
	ElemGroup elG;
	
	Coupling = new vector<intVec*>(s->NSubsteps); 
	for(int i=0; i<s->NSubsteps; i++) Coupling->at(i) = new intVec(n);

	for(int k = 0; k< s->NSubsteps; k++){		
		
		for(int e=0; e< d->Nelem; e++){	
			int iGroup = (*d->eType)(e);
			
			elG = s->substeps->at(k).elGroup->at(iGroup);
			
			(*CouplingE)=0;
			
			if(elG.FlagElemLib<0){
				getLocalSymbolic(&elG.ElemLib, CouplingE->v, elG.CommonPar->v,&d->Ndof, &d->Ndim, &MaxLRows);
			}
			
			assemblyCoupling(Coupling->at(k),CouplingE,e);
		}
		
		for(int i=0; i<n; i++) (*Coupling->at(k))(i) = (*Coupling->at(k))(i) - 1; 
	}
}

void Solver:: linsolve(){
	solve(*Sol,*A,*b);
}
	 
void Solver :: enforceDirichlet(SGPrealVec *vDirich, intVec *flagDirich, double pen = 99999.9){

	for(int i = 0 ; i<b->size(); i++){
		if((*flagDirich)(i)==-1){
			(*b)(i) = pen*(*vDirich)(i);
			(*A)(i,i) = pen*1.0;
			for(int j = 0 ; j<b->size(); j++){ if(j!=i) (*A)(i,j) = 0.0;}
		}
	}
}
//~ 
void Solver :: assembly(SGPrealMat *A_L, SGPrealVec *b_L, int e){
	
	int ip,jp,kip,kjp;
	
	int iShiftElem = (*d->eNNod_acc)(e);
	
	for(int i=0; i<(*d->eNNod)(e); i++){ for(int ii=0; ii<d->Ndof; ii++){
		ip = i*d->Ndof + ii;
		kip = (*d->Elem)(iShiftElem + i)*d->Ndof + ii;
		(*b)(kip) += (*b_L)(ip);
		for(int j=0; j<(*d->eNNod)(e); j++){ for(int jj=0; jj<d->Ndof; jj++){
			jp = j*d->Ndof + jj;
			kjp = (*d->Elem)(iShiftElem + j)*d->Ndof + jj;
			(*A)(kip,kjp) += (*A_L)(ip,jp);			
		}
		}
	}
	}
}

void Solver :: assemblyCoupling(intVec *Cvec, intMat *Cmat_L, int e){
	
	int ip,kip;
	int iShiftElem = (*d->eNNod_acc)(e);
	
	for(int i=0; i<(*d->eNNod)(e); i++){ for(int ii=0; ii<d->Ndof; ii++){
		ip = i*d->Ndof + ii;
		kip = (*d->Elem)(iShiftElem + i)*d->Ndof + ii;
		if((*Cmat_L)(ip,ip) == 1){
			(*Cvec)(kip) = 1;
		}
	}
	}
}

void Solver :: globalToLocal(SGPrealVec *XLL,SGPrealVec *Sol0E, SGPrealVec *Sol1E, SGPrealVec *ParamE, int e){
		int ipDim,kpDim,ipDof,kpDof;
		int iShiftElem = (*d->eNNod_acc)(e);
		int iShiftMat = (*d->eMat)(e);
		int iShiftParam = (*d->eParamSize_acc)(iShiftMat);
	
		(*XLL) = 0.0;
		(*Sol0E) = 0.0;
		(*Sol1E) = 0.0;
		(*ParamE) = 0.0;
		 
		for(int i=0;i<(*d->eNNod)(e);i++){
			ipDim = i*d->Ndim;
			kpDim = (*d->Elem)(iShiftElem + i)*d->Ndim;
			ipDof = i*d->Ndof;
			kpDof = (*d->Elem)(iShiftElem + i)*d->Ndof;
					
			for(int j=0 ; j<d->Ndim; j++) (*XLL)(ipDim + j) = (*d->X)(kpDim +j);
			
			for(int j=0 ; j<d->Ndof; j++){
				(*Sol0E)(ipDof + j) = (*d->Sol0)(kpDof + j);
				(*Sol1E)(ipDof + j) = (*d->Sol1)(kpDof + j);
			}
		}
		
		for(int k=0; k<(*d->eParamSize)(iShiftMat) ; k++) (*ParamE)(k) = (*d->Param)(iShiftParam + k);
		
}

void Solver :: localToGlobal(SGPrealVec *ParamE, int e){
		int iShiftMat = (*d->eMat)(e);
		int iShiftParam = (*d->eParamSize_acc)(iShiftMat);
			 	
		for(int k=0; k<(*d->eParamSize)(iShiftMat) ; k++) (*d->Param)(iShiftParam + k) = (*ParamE)(k);
}

void Solver :: run(){
	int timeStep, iGroup;
	int id_Elem_Family, MaxLRows, iDofT, NodElt;
	intVec *JParam;
	double DTm, Time;
	SGPrealVec *BE, *XLL, *Sol0E, *Sol1E, *CommonParE, *ParamE;
	SGPrealMat *AE; 
	ElemGroup elG;
	double error, tol, emax = 99999.0;
	int itNL, maxit;
	bool isNL;
	
	NodElt = d->eNNod->max();
	MaxLRows = NodElt*d->Ndof;

	AE = new SGPrealMat(MaxLRows,MaxLRows);
	BE = new SGPrealVec(MaxLRows);
	XLL = new SGPrealVec(NodElt*d->Ndim);
	Sol0E = new SGPrealVec(MaxLRows);
	Sol1E = new SGPrealVec(MaxLRows);
	ParamE = new SGPrealVec(d->eParamSize->max()); 
	JParam = new intVec(3); // temporary
	
	Time = s->Tini; 
	DTm = s->DelT; 	
	timeStep = 0;
	
	d->writeSGPdataout(timeStep,Time,DTm);
	
	while(Time < s->Tmax){
		DEBUG{cout << "=> Beginning Timestep = " << timeStep << ", time = " << Time << "  =====\n " << endl;}
		
		for(int k = 0; k< s->NSubsteps; k++){
			
			error = emax;
			itNL = 0;
			isNL = s->substeps->at(k).isNL;
			if(isNL){
				maxit = s->substeps->at(k).maxit;
				tol = s->substeps->at(k).tol;
			}
			else{
				maxit = 1; // force to stop with just one iteration
				tol = 0.0;
			}
			
			DEBUG{cout << "===> Beginning SubStep = " << k << " ===========" << endl;}
			
			DEBUG{if(isNL) cout << "===> Starting NonLinear Convergence, tol = " << tol << endl;}
			while(error>tol && itNL<maxit){
				
				DEBUG{if(isNL) cout << "===> Starting NonLinear Internal Loop, it = " << itNL << " =======" << endl;} 
			
				A->fill(0.0);
				b->fill(0.0);
					
				//~ cout << "=====> Assemblying Global Matrices ================ " << endl;
				for(int e=0; e< d->Nelem; e++){
					
					globalToLocal(XLL,Sol0E,Sol1E,ParamE,e);
					
					iGroup = (*d->eType)(e);
					
					elG = s->substeps->at(k).elGroup->at(iGroup);
					
					(*AE)=0.0;
					(*BE)=0.0;
					
					if(elG.FlagElemLib<0){
						getLocalMatrix(&elG.ElemLib, AE->v, BE->v, &MaxLRows, XLL->v, &d->Ndim, &d->Ndof, &NodElt, Sol0E->v, Sol1E->v, 
										elG.CommonPar->v, ParamE->v, JParam->v, &s->DelT, &DTm, &Time);
					}
					
					assembly(AE,BE,e);
					
					localToGlobal(ParamE,e);
				}
				
				DEBUG{cout << "=====> Solving Linear System ================ " << endl;}
				enforceDirichlet( d->vDirich->at(k), d->flagDirich->at(k) , 1.0);
				enforceDirichlet( d->Sol1, Coupling->at(k) , 1.0);
					
				linsolve();
				if(isNL){
					error = computeError(*Sol,*d->Sol1);
					DEBUG{cout << "====> Ending NonLinear Internal Loop , it= " << itNL << " , error= " << error << " =========="<< endl;}	
				} 
				
				itNL ++;
				*d->Sol1 = *Sol; // uses the operator overloading
			}
			DEBUG{if(isNL) cout << "====> Convergence Achieved with it= " << itNL-1 << " , error= " << error << " =========="<< endl;}
			
			DEBUG{cout << "===> Ending SubStep = " << k << " =============== \n " << endl;}
		}
		
		DEBUG{cout << "=> Ending Timestep = " << timeStep << ", time = " << Time << " ===== \n" << endl;}
		
		Time += s->DelT;
		timeStep++;
		
		*d->Sol0 = *d->Sol1; // uses the operator overloading
		
		DEBUG{cout << "=> Writing Results in File ====== \n\n" << endl;}
		d->writeSGPdataout(timeStep,Time,DTm);
		d->writeSGPparam();
	}
}



#endif
