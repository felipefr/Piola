# file CMakeLists.txt
 
# esta primeira linha não é obrigatória. Mas o CMake às vezes reclama da ausência dela
cmake_minimum_required(VERSION 2.6)
 
# definição do nome do project. Não precisa ser obrigatoriamente o nome do diretório
project(Piola)

SET(LIBMATERIALS /home/felipefr/github/Piola/fortranInterface/include)
SET(LIBMATERIALS_LIB /home/felipefr/github/Piola/fortranInterface/libroutines.so)
SET(ARMADILLO_DIR /home/felipefr/.programFiles/armadillo/armadillo-6.500.5/include)
SET(ARMADILLO_LIB_SO /home/felipefr/.programFiles/armadillo/armadillo-6.500.5/libarmadillo.so)

SET(Directoires ../src/)
SET(fileExtensions .c .cpp .cc .h .hpp) 
SET(src)
foreach(d ${Directoires})
	foreach(e ${fileExtensions})
		file(GLOB dSources ${d}/*${e})
		LIST(APPEND src ${dSources})
	endforeach(e)
endforeach(d)

include_directories(${LIBMATERIALS} ${ARMADILLO_DIR})

add_executable(Piola.x ${src})

target_link_libraries(Piola.x ${LIBMATERIALS_LIB} ${ARMADILLO_LIB} lapack) 
