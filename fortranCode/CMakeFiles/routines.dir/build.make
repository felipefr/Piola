# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/felipefr/.programFiles/cmake-3.20.3-linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/felipefr/.programFiles/cmake-3.20.3-linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/felipefr/github/Piola/fortranCode

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/felipefr/github/Piola/fortranCode

# Include any dependencies generated for this target.
include CMakeFiles/routines.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/routines.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/routines.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/routines.dir/flags.make

CMakeFiles/routines.dir/routines.f90.o: CMakeFiles/routines.dir/flags.make
CMakeFiles/routines.dir/routines.f90.o: routines.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felipefr/github/Piola/fortranCode/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/routines.dir/routines.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/felipefr/github/Piola/fortranCode/routines.f90 -o CMakeFiles/routines.dir/routines.f90.o

CMakeFiles/routines.dir/routines.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/routines.dir/routines.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/felipefr/github/Piola/fortranCode/routines.f90 > CMakeFiles/routines.dir/routines.f90.i

CMakeFiles/routines.dir/routines.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/routines.dir/routines.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/felipefr/github/Piola/fortranCode/routines.f90 -o CMakeFiles/routines.dir/routines.f90.s

# Object files for target routines
routines_OBJECTS = \
"CMakeFiles/routines.dir/routines.f90.o"

# External object files for target routines
routines_EXTERNAL_OBJECTS =

libroutines.so: CMakeFiles/routines.dir/routines.f90.o
libroutines.so: CMakeFiles/routines.dir/build.make
libroutines.so: /home/felipefr/github/Piola/gpMaterials_beta/libgpmaterials.so
libroutines.so: CMakeFiles/routines.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/felipefr/github/Piola/fortranCode/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran shared library libroutines.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/routines.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/routines.dir/build: libroutines.so
.PHONY : CMakeFiles/routines.dir/build

CMakeFiles/routines.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/routines.dir/cmake_clean.cmake
.PHONY : CMakeFiles/routines.dir/clean

CMakeFiles/routines.dir/depend:
	cd /home/felipefr/github/Piola/fortranCode && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/felipefr/github/Piola/fortranCode /home/felipefr/github/Piola/fortranCode /home/felipefr/github/Piola/fortranCode /home/felipefr/github/Piola/fortranCode /home/felipefr/github/Piola/fortranCode/CMakeFiles/routines.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/routines.dir/depend

