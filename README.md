# Piola
Piola aims at mimetising SolverGP to provide a simple, C++, opensource alternative to it. I developed during my PhD at LNCC (2015-2019) for debug purposes. Please consider compilating a version of GPMaterials in a .so file, to be able to call the library from Piola. This version specifically relies on a valid Armadillo instalation for linearalgebra routines.


i) sudo apt-get install libarmadillo-dev
ii) in the home directory of Piola, access the fortranCode directory
iii) change path to a valid library (static or shared) of formulations according to the interface (e.g. as the https://github.com/felipefr/gpmaterials)
iv) cmake .
v) make
vi) cd ../build
vii) cmake .
viii) make

TODO: some errors arising with lapack version armadillo
