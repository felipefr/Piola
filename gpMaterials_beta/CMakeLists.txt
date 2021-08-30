cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)
 
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(gpmaterials)

ADD_DEFINITIONS(-fPIC -g -fbacktrace)

SET(Directoires src/tools src/otherMaterials src/damage src/Multiscale src/FiniteStrain src/fibres src/fe2 src/identification src/EFEM)
SET(fileExtensions .f .for .f90 .F .F90 .FOR) 
SET(src)
foreach(d ${Directoires})
	foreach(e ${fileExtensions})
		file(GLOB dSources ${d}/*${e})
		LIST(APPEND src ${dSources})
	endforeach(e)
endforeach(d)

add_library(gpmaterials SHARED ${src})
