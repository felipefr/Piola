# Piola
Piola aims at mimetising SolverGP to provide a simple, C++, opensource alternative to it. I developed during my PhD at LNCC (2015-2019) for debug purposes. Please consider compilating a version of GPMaterials in a .so file, to be able to call the library from Piola. This version specifically relies on a valid Armadillo instalation for linearalgebra routines.


i) Install armadillo version 6.500 with the necessary dependencies (maybe other versions also work)

ii) cd myelements/ (in case of having your own library of elements jump to step iv), e.g the one found in https://github.com/felipefr/gpmaterials), but don't forget in changing CMakeLists.txt in the fortranInterface folder)

iii) cmake . ; make

iv) cd ../fortranInterface 

v) cmake . ; make

vi) cd ../build

vii) cmake .; make

