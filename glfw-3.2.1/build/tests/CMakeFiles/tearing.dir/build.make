# CMAKE generated file: DO NOT EDIT!
# Generated by "MSYS Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = "/C/Program Files/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "/C/Program Files/CMake/bin/cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/tearing.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/tearing.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/tearing.dir/flags.make

tests/CMakeFiles/tearing.dir/tearing.c.obj: tests/CMakeFiles/tearing.dir/flags.make
tests/CMakeFiles/tearing.dir/tearing.c.obj: tests/CMakeFiles/tearing.dir/includes_C.rsp
tests/CMakeFiles/tearing.dir/tearing.c.obj: ../tests/tearing.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/tearing.dir/tearing.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tearing.dir/tearing.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/tearing.c

tests/CMakeFiles/tearing.dir/tearing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/tearing.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/tearing.c > CMakeFiles/tearing.dir/tearing.c.i

tests/CMakeFiles/tearing.dir/tearing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/tearing.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/tearing.c -o CMakeFiles/tearing.dir/tearing.c.s

tests/CMakeFiles/tearing.dir/__/deps/getopt.c.obj: tests/CMakeFiles/tearing.dir/flags.make
tests/CMakeFiles/tearing.dir/__/deps/getopt.c.obj: tests/CMakeFiles/tearing.dir/includes_C.rsp
tests/CMakeFiles/tearing.dir/__/deps/getopt.c.obj: ../deps/getopt.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object tests/CMakeFiles/tearing.dir/__/deps/getopt.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tearing.dir/__/deps/getopt.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c

tests/CMakeFiles/tearing.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/__/deps/getopt.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c > CMakeFiles/tearing.dir/__/deps/getopt.c.i

tests/CMakeFiles/tearing.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/__/deps/getopt.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c -o CMakeFiles/tearing.dir/__/deps/getopt.c.s

tests/CMakeFiles/tearing.dir/__/deps/glad.c.obj: tests/CMakeFiles/tearing.dir/flags.make
tests/CMakeFiles/tearing.dir/__/deps/glad.c.obj: tests/CMakeFiles/tearing.dir/includes_C.rsp
tests/CMakeFiles/tearing.dir/__/deps/glad.c.obj: ../deps/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object tests/CMakeFiles/tearing.dir/__/deps/glad.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tearing.dir/__/deps/glad.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c

tests/CMakeFiles/tearing.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/__/deps/glad.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c > CMakeFiles/tearing.dir/__/deps/glad.c.i

tests/CMakeFiles/tearing.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/__/deps/glad.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c -o CMakeFiles/tearing.dir/__/deps/glad.c.s

# Object files for target tearing
tearing_OBJECTS = \
"CMakeFiles/tearing.dir/tearing.c.obj" \
"CMakeFiles/tearing.dir/__/deps/getopt.c.obj" \
"CMakeFiles/tearing.dir/__/deps/glad.c.obj"

# External object files for target tearing
tearing_EXTERNAL_OBJECTS =

tests/tearing.exe: tests/CMakeFiles/tearing.dir/tearing.c.obj
tests/tearing.exe: tests/CMakeFiles/tearing.dir/__/deps/getopt.c.obj
tests/tearing.exe: tests/CMakeFiles/tearing.dir/__/deps/glad.c.obj
tests/tearing.exe: tests/CMakeFiles/tearing.dir/build.make
tests/tearing.exe: src/libglfw3dll.a
tests/tearing.exe: tests/CMakeFiles/tearing.dir/linklibs.rsp
tests/tearing.exe: tests/CMakeFiles/tearing.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable tearing.exe"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && "/C/Program Files/CMake/bin/cmake.exe" -E remove -f CMakeFiles/tearing.dir/objects.a
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/ar.exe cr CMakeFiles/tearing.dir/objects.a @CMakeFiles/tearing.dir/objects1.rsp
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe   -mwindows -Wl,--whole-archive CMakeFiles/tearing.dir/objects.a -Wl,--no-whole-archive  -o tearing.exe -Wl,--out-implib,libtearing.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/tearing.dir/linklibs.rsp

# Rule to build all files generated by this target.
tests/CMakeFiles/tearing.dir/build: tests/tearing.exe

.PHONY : tests/CMakeFiles/tearing.dir/build

tests/CMakeFiles/tearing.dir/clean:
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/tearing.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/tearing.dir/clean

tests/CMakeFiles/tearing.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MSYS Makefiles" /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1 /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests/CMakeFiles/tearing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/tearing.dir/depend

