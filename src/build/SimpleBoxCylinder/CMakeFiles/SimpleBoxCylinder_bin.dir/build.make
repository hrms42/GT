# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/HRMS/Desktop/GT/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/HRMS/Desktop/GT/src/build

# Include any dependencies generated for this target.
include SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/depend.make

# Include the progress variables for this target.
include SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/progress.make

# Include the compile flags for this target's objects.
include SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/flags.make

SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o: SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/flags.make
SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o: ../SimpleBoxCylinder/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/HRMS/Desktop/GT/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o"
	cd /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o -c /Users/HRMS/Desktop/GT/src/SimpleBoxCylinder/main.cpp

SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.i"
	cd /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/HRMS/Desktop/GT/src/SimpleBoxCylinder/main.cpp > CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.i

SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.s"
	cd /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/HRMS/Desktop/GT/src/SimpleBoxCylinder/main.cpp -o CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.s

# Object files for target SimpleBoxCylinder_bin
SimpleBoxCylinder_bin_OBJECTS = \
"CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o"

# External object files for target SimpleBoxCylinder_bin
SimpleBoxCylinder_bin_EXTERNAL_OBJECTS =

SimpleBoxCylinder_bin: SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/main.cpp.o
SimpleBoxCylinder_bin: SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/build.make
SimpleBoxCylinder_bin: glad/libglad.a
SimpleBoxCylinder_bin: glfw/src/libglfw3.a
SimpleBoxCylinder_bin: SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/HRMS/Desktop/GT/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../SimpleBoxCylinder_bin"
	cd /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SimpleBoxCylinder_bin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/build: SimpleBoxCylinder_bin

.PHONY : SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/build

SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/clean:
	cd /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder && $(CMAKE_COMMAND) -P CMakeFiles/SimpleBoxCylinder_bin.dir/cmake_clean.cmake
.PHONY : SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/clean

SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/depend:
	cd /Users/HRMS/Desktop/GT/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/HRMS/Desktop/GT/src /Users/HRMS/Desktop/GT/src/SimpleBoxCylinder /Users/HRMS/Desktop/GT/src/build /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder /Users/HRMS/Desktop/GT/src/build/SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : SimpleBoxCylinder/CMakeFiles/SimpleBoxCylinder_bin.dir/depend

