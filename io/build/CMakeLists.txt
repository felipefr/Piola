# esta primeira linha não é obrigatória. Mas o CMake às vezes reclama da ausência dela
cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)
  
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(converter2018)

ADD_DEFINITIONS(-fPIC -g -fbacktrace)

#~ # OPTION 2 = no library making
SET(Directoires ../src)
SET(fileExtensions .f .for .f90 .F .F90 .FOR) 
SET(src)
foreach(d ${Directoires})
	foreach(e ${fileExtensions})
		file(GLOB dSources ${d}/*${e})
		LIST(APPEND src ${dSources})
	endforeach(e)
endforeach(d)

add_executable(converter.x ${src})

