# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.0.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.0.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/takashi/cbrc/ractip_hom/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/takashi/cbrc/ractip_hom/src

# Include any dependencies generated for this target.
include centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/depend.make

# Include the progress variables for this target.
include centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/progress.make

# Include the compile flags for this target's objects.
include centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/flags.make

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/flags.make
centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o: centroidalifold/probconsRNA/wrapper.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/takashi/cbrc/ractip_hom/src/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Probcons.dir/wrapper.cpp.o -c /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA/wrapper.cpp

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Probcons.dir/wrapper.cpp.i"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA/wrapper.cpp > CMakeFiles/Probcons.dir/wrapper.cpp.i

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Probcons.dir/wrapper.cpp.s"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA/wrapper.cpp -o CMakeFiles/Probcons.dir/wrapper.cpp.s

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.requires:
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.requires

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.provides: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.requires
	$(MAKE) -f centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/build.make centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.provides.build
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.provides

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.provides.build: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o

# Object files for target Probcons
Probcons_OBJECTS = \
"CMakeFiles/Probcons.dir/wrapper.cpp.o"

# External object files for target Probcons
Probcons_EXTERNAL_OBJECTS =

centroidalifold/probconsRNA/libProbcons.a: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o
centroidalifold/probconsRNA/libProbcons.a: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/build.make
centroidalifold/probconsRNA/libProbcons.a: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libProbcons.a"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && $(CMAKE_COMMAND) -P CMakeFiles/Probcons.dir/cmake_clean_target.cmake
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Probcons.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/build: centroidalifold/probconsRNA/libProbcons.a
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/build

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/requires: centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/wrapper.cpp.o.requires
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/requires

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/clean:
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA && $(CMAKE_COMMAND) -P CMakeFiles/Probcons.dir/cmake_clean.cmake
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/clean

centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/depend:
	cd /Users/takashi/cbrc/ractip_hom/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA /Users/takashi/cbrc/ractip_hom/src/centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : centroidalifold/probconsRNA/CMakeFiles/Probcons.dir/depend

