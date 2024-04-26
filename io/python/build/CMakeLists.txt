cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)


# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(converter)

ADD_DEFINITIONS(-fPIC)

SET(src homogenisationPosProcessing.f90 utils.f90)

add_library(converter STATIC ${src})
