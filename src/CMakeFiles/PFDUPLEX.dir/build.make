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
include CMakeFiles/PFDUPLEX.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PFDUPLEX.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PFDUPLEX.dir/flags.make

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o: CMakeFiles/PFDUPLEX.dir/flags.make
CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o: pf_duplex.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/takashi/cbrc/ractip_hom/src/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o   -c /Users/takashi/cbrc/ractip_hom/src/pf_duplex.c

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/PFDUPLEX.dir/pf_duplex.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/takashi/cbrc/ractip_hom/src/pf_duplex.c > CMakeFiles/PFDUPLEX.dir/pf_duplex.c.i

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/PFDUPLEX.dir/pf_duplex.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/takashi/cbrc/ractip_hom/src/pf_duplex.c -o CMakeFiles/PFDUPLEX.dir/pf_duplex.c.s

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.requires:
.PHONY : CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.requires

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.provides: CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.requires
	$(MAKE) -f CMakeFiles/PFDUPLEX.dir/build.make CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.provides.build
.PHONY : CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.provides

CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.provides.build: CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o

# Object files for target PFDUPLEX
PFDUPLEX_OBJECTS = \
"CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o"

# External object files for target PFDUPLEX
PFDUPLEX_EXTERNAL_OBJECTS =

libPFDUPLEX.a: CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o
libPFDUPLEX.a: CMakeFiles/PFDUPLEX.dir/build.make
libPFDUPLEX.a: CMakeFiles/PFDUPLEX.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libPFDUPLEX.a"
	$(CMAKE_COMMAND) -P CMakeFiles/PFDUPLEX.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PFDUPLEX.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PFDUPLEX.dir/build: libPFDUPLEX.a
.PHONY : CMakeFiles/PFDUPLEX.dir/build

CMakeFiles/PFDUPLEX.dir/requires: CMakeFiles/PFDUPLEX.dir/pf_duplex.c.o.requires
.PHONY : CMakeFiles/PFDUPLEX.dir/requires

CMakeFiles/PFDUPLEX.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PFDUPLEX.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PFDUPLEX.dir/clean

CMakeFiles/PFDUPLEX.dir/depend:
	cd /Users/takashi/cbrc/ractip_hom/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src /Users/takashi/cbrc/ractip_hom/src/CMakeFiles/PFDUPLEX.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PFDUPLEX.dir/depend
