// implements the SimulSettings class - Felipe Figueredo Rocha 
#ifndef _simulSettings_hpp
#define _simulSettings_hpp

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

using namespace std;

class ElemGroup{
	public:
	int ElemLib, FlagElemLib;
	SGPrealVec *CommonPar;
	ElemGroup() {}
};

class SubStep{
	public:
	bool isNL;
	double tol;
	int maxit;
	vector<ElemGroup> *elGroup; 
	SubStep() {}
};

class SimulSettings{
	public:
	// see constructor for further explanation about the variabls
	int NSubsteps, MaxElemLib,MaxLCommonPar,MaxNodEl; // sizes
	double DelT, Tini, Tmax;
	vector<SubStep> *substeps;
	SimulSettings() {}
	void loadSettings(); // read a mesh  
};

void SimulSettings :: loadSettings(){
	
	ifstream file("Basparam.txt") ;
	string line;
	
	findString(file , "*TimeStep"); 
	file >> DelT >> Tini >> Tmax ;

	findString(file , "*ElementLibraryControl"); 
	file >> NSubsteps >> MaxElemLib >> MaxLCommonPar >> MaxNodEl; 
	
	substeps = new vector<SubStep>(NSubsteps);
	
	for(int i=0; i<NSubsteps; i++){
	
		substeps->at(i).elGroup = new vector<ElemGroup>(MaxElemLib);
		
		findString(file , "*SubStep"); 	
		file >> substeps->at(i).isNL;
		if(substeps->at(i).isNL) file >> substeps->at(i).tol >> substeps->at(i).maxit ;
		 
		for(int j=0; j<MaxElemLib; j++){ 
			
			file >> substeps->at(i).elGroup->at(j).ElemLib >> substeps->at(i).elGroup->at(j).FlagElemLib;
			substeps->at(i).elGroup->at(j).CommonPar = new SGPrealVec(MaxLCommonPar);
		}
	
		for(int j=0; j<MaxElemLib; j++) substeps->at(i).elGroup->at(j).CommonPar->readNumberBlock(file);
	}
	
	file.close();
}	


#endif
