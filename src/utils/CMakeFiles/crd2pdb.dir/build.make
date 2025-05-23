# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yamauchi/cafemol

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yamauchi/cafemol

# Include any dependencies generated for this target.
include src/utils/CMakeFiles/crd2pdb.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/utils/CMakeFiles/crd2pdb.dir/compiler_depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/crd2pdb.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/crd2pdb.dir/flags.make

src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.o: src/utils/CMakeFiles/crd2pdb.dir/flags.make
src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.o: src/utils/crd2pdb.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yamauchi/cafemol/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.o"
	cd /home/yamauchi/cafemol/src/utils && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/yamauchi/cafemol/src/utils/crd2pdb.F90 -o CMakeFiles/crd2pdb.dir/crd2pdb.F90.o

src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/crd2pdb.dir/crd2pdb.F90.i"
	cd /home/yamauchi/cafemol/src/utils && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/yamauchi/cafemol/src/utils/crd2pdb.F90 > CMakeFiles/crd2pdb.dir/crd2pdb.F90.i

src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/crd2pdb.dir/crd2pdb.F90.s"
	cd /home/yamauchi/cafemol/src/utils && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/yamauchi/cafemol/src/utils/crd2pdb.F90 -o CMakeFiles/crd2pdb.dir/crd2pdb.F90.s

# Object files for target crd2pdb
crd2pdb_OBJECTS = \
"CMakeFiles/crd2pdb.dir/crd2pdb.F90.o"

# External object files for target crd2pdb
crd2pdb_EXTERNAL_OBJECTS =

bin/crd2pdb: src/utils/CMakeFiles/crd2pdb.dir/crd2pdb.F90.o
bin/crd2pdb: src/utils/CMakeFiles/crd2pdb.dir/build.make
bin/crd2pdb: src/libcafemol.a
bin/crd2pdb: src/utils/CMakeFiles/crd2pdb.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yamauchi/cafemol/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable ../../bin/crd2pdb"
	cd /home/yamauchi/cafemol/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/crd2pdb.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/utils/CMakeFiles/crd2pdb.dir/build: bin/crd2pdb
.PHONY : src/utils/CMakeFiles/crd2pdb.dir/build

src/utils/CMakeFiles/crd2pdb.dir/clean:
	cd /home/yamauchi/cafemol/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/crd2pdb.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/crd2pdb.dir/clean

src/utils/CMakeFiles/crd2pdb.dir/depend:
	cd /home/yamauchi/cafemol && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yamauchi/cafemol /home/yamauchi/cafemol/src/utils /home/yamauchi/cafemol /home/yamauchi/cafemol/src/utils /home/yamauchi/cafemol/src/utils/CMakeFiles/crd2pdb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/crd2pdb.dir/depend

