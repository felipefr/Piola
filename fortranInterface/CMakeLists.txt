# esta primeira linha não é obrigatória. Mas o CMake às vezes reclama da ausência dela
cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)
  
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(routines)

# SET(LIBMATERIALS_LIB /home/felipefr/github/gpmaterials/build/libgpmaterials.so)
SET(LIBMATERIALS_LIB /home/felipefr/github/Piola/myelements/libmyelements.so)

ADD_DEFINITIONS(-fPIC -g -fbacktrace)

SET(src routines.f90)

add_library(routines SHARED ${src})

target_link_libraries(routines ${LIBMATERIALS_LIB})
