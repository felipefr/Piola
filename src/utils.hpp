// some utilitary routines - Felipe Figueredo Rocha
#ifndef _utils_hpp
#define _utils_hpp

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>  
#include <sstream>
#include <cmath>
#include "linAlg.hpp"

using namespace std;

const SGPreal PI=4.0*atan(1.0);

// stops when finds a the string searched
void findString(ifstream &file, string search_str) 
{
	string line ;
        
	while( getline( file, line ) ){
		if( line.find(search_str) != string::npos ) return;
	}
}

#endif
