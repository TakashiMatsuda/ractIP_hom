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
include centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/depend.make

# Include the progress variables for this target.
include centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/progress.make

# Include the compile flags for this target's objects.
include centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/flags.make

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/flags.make
centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o: centroidalifold/folding_engine.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/takashi/cbrc/ractip_hom/src/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o -c /Users/takashi/cbrc/ractip_hom/src/centroidalifold/folding_engine.cpp

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.i"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/takashi/cbrc/ractip_hom/src/centroidalifold/folding_engine.cpp > CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.i

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.s"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/takashi/cbrc/ractip_hom/src/centroidalifold/folding_engine.cpp -o CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.s

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.requires:
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.requires

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.provides: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.requires
	$(MAKE) -f centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/build.make centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.provides.build
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.provides

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.provides.build: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o

# Object files for target FOLDING_ENGINE
FOLDING_ENGINE_OBJECTS = \
"CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o"

# External object files for target FOLDING_ENGINE
FOLDING_ENGINE_EXTERNAL_OBJECTS =

centroidalifold/libFOLDING_ENGINE.a: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o
centroidalifold/libFOLDING_ENGINE.a: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/build.make
centroidalifold/libFOLDING_ENGINE.a: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libFOLDING_ENGINE.a"
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && $(CMAKE_COMMAND) -P CMakeFiles/FOLDING_ENGINE.dir/cmake_clean_target.cmake
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FOLDING_ENGINE.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/build: centroidalifold/libFOLDING_ENGINE.a
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/build

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/requires: centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/folding_engine.cpp.o.requires
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/requires

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/clean:
	cd /Users/takashi/cbrc/ractip_hom/src/centroidalifold && $(CMAKE_COMMAND) -P CMakeFiles/FOLDING_ENGINE.dir/cmake_clean.cmake
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/clean

centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/depend:
	cd /Users/takashi/cbrc/ractip_hom/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src/centroidalifold /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src/centroidalifold /Users/takashi/cbrc/ractip_hom/src/centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : centroidalifold/CMakeFiles/FOLDING_ENGINE.dir/depend
