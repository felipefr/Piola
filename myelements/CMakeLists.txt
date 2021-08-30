cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)
 
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(myelements)

ADD_DEFINITIONS(-fPIC -g -fbacktrace)

SET(Directoires src/tools src/otherMaterials src/FiniteStrain)
SET(fileExtensions .f .for .f90 .F .F90 .FOR) 
SET(src)
foreach(d ${Directoires})
	foreach(e ${fileExtensions})
		file(GLOB dSources ${d}/*${e})
		LIST(APPEND src ${dSources})
	endforeach(e)
endforeach(d)

add_library(myelements SHARED ${src})
