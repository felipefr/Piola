// implements the Data class - Felipe Figueredo Rocha 
#ifndef _data_hpp
#define _data_hpp

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

class Data{
	public:
	// see constructor for further explanation about the variabls
	int Nelem, Nnodes, Ndim, Ndof; // sizes
	intVec *Elem, *eType, *eMat, *eNNod, *eNNod_acc, *eParamSize, *eParamSize_acc; 
	vector<intVec*> *flagDirich; 
	vector<SGPrealVec*> *vDirich;
	SGPrealVec *X, *Sol0,*Sol1,*Param;
	Data() {}
	void loadSGPmesh(); // read a mesh  
	void loadSGPdirich(int nSubSteps); // read a Dirichlet
	void loadSGPinifile(); // read a Inifile  
	void loadSGPparam(); // read a Param
	void writeSGPdataout(int timeStep,double time, double dt); // write results
	void writeSGPparam(); // write results
};

void Data:: writeSGPdataout(int timeStep,double time, double dt){
	ofstream file("DataOut.txt",ios_base::app) ;
	
	file << "*Time" << endl;
	file << timeStep << " " << time << " " << dt << endl;
	Sol1->writeNumberBlockRect(file,Ndof);
	
	file.close();
}

void Data :: loadSGPmesh(){
	int dummy;
	int NGroups;
	intVec *NeGroup;

	ifstream file("Mesh.txt") ;
	string line;
	
	findString(file , "*NODAL DOFs"); 
	file >> Ndof;
	findString(file , "*DIMEN"); 
	file >> Ndim;
	
	findString(file , "*COORDINATES" );
	file >> Nnodes;
	X = new SGPrealVec(Ndim*Nnodes);
	X->readNumberBlock(file);
	
	findString(file , "*ELEMENT GROUPS");
	file >> NGroups;
	NeGroup = new intVec(NGroups);
	for(int i=0; i<NGroups; i++){ file >> dummy; file >> NeGroup->v[i]; getline(file,line); } 
	 
	Nelem = NeGroup->sum();
	
	eNNod = new intVec(Nelem);
	eNNod->readNumberBlock(file);
	
	eNNod_acc = new intVec(Nelem);
	for(int i=1; i< Nelem ; i++) (*eNNod_acc)(i) = (*eNNod_acc)(i-1) + (*eNNod)(i-1);
	
	Elem = new intVec(eNNod->sum());
	findString(file , "*INCIDENCE");
	Elem->readNumberBlock(file);
	for(int i = 0; i<Elem->n ; i++) (*Elem)(i) -= 1;
	
	eType = new intVec(Nelem);
	findString(file , "*ELEMENT TYPE");
	eType->readNumberBlock(file);
	for(int i = 0; i<Nelem ; i++) (*eType)(i) -= 1;
		
	eMat = new intVec(Nelem);
	findString(file , "*ELEMENT MAT");
	eMat->readNumberBlock(file);
	for(int i = 0; i<Nelem ; i++) (*eMat)(i) -= 1;
	
	file.close();
}	

void Data :: loadSGPdirich(int nSubSteps){
	int dummy;

	ifstream file("Mesh.txt") ;
	string line;
	
	flagDirich = new vector<intVec*>(nSubSteps);
	vDirich = new vector<SGPrealVec*>(nSubSteps);
	
	
	for(int i=0; i<nSubSteps; i++){
		flagDirich->at(i) = new intVec(Nnodes*Ndof);
		vDirich->at(i) = new SGPrealVec(Nnodes*Ndof);
	}
		
	findString(file , "*DIRICHLET CONDITIONS");
	
	for(int i=0; i<nSubSteps; i++) flagDirich->at(i)->readNumberBlock(file);
	for(int i=0; i<nSubSteps; i++) vDirich->at(i)->readNumberBlock(file);
	
	file.close();
}	

void Data :: loadSGPinifile(){

	ifstream file("IniFile.txt");
	
	Sol0 = new SGPrealVec(Nnodes*Ndof);
	Sol1 = new SGPrealVec(Nnodes*Ndof);
	findString(file , "*Initial Conditions");
	Sol0->readNumberBlock(file);
	
	file.close();
}

void Data :: loadSGPparam(){
	int NparamGroups, totalParamSize;

	ifstream file("Param.txt");
	
	findString(file , "*Parameter Groups");
	file >> NparamGroups;
	
	eParamSize = new intVec(NparamGroups);
	
	findString(file , "*Real Parameters");
	
	eParamSize->readNumberBlock(file);
	
	eParamSize_acc = new intVec(NparamGroups);
	for(int i=1; i< NparamGroups ; i++) (*eParamSize_acc)(i) = (*eParamSize_acc)(i-1) + (*eParamSize)(i-1);
	
	totalParamSize = eParamSize->sum();
	
	Param = new SGPrealVec(totalParamSize);
	Param->readNumberBlock(file);
	
	file.close();
}		

void Data:: writeSGPparam(){
	ofstream file("Param.txt",ios_base::app) ;
	
	file << "*Parameter Groups" << endl;
	file << eParamSize->n  << endl;
	
	file << "*Real Parameters" << endl;
	eParamSize->writeNumberBlock(file);
	Param->writeNumberBlock(file);
	
	file << "*Integer Parameters" << endl;
	for(int i=0; i< eParamSize->n; i++) file << 0 << endl;

	file.close();
}
#endif
