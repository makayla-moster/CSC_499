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
include tests/CMakeFiles/msaa.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/msaa.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/msaa.dir/flags.make

tests/CMakeFiles/msaa.dir/msaa.c.obj: tests/CMakeFiles/msaa.dir/flags.make
tests/CMakeFiles/msaa.dir/msaa.c.obj: tests/CMakeFiles/msaa.dir/includes_C.rsp
tests/CMakeFiles/msaa.dir/msaa.c.obj: ../tests/msaa.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object tests/CMakeFiles/msaa.dir/msaa.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/msaa.dir/msaa.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/msaa.c

tests/CMakeFiles/msaa.dir/msaa.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/msaa.dir/msaa.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/msaa.c > CMakeFiles/msaa.dir/msaa.c.i

tests/CMakeFiles/msaa.dir/msaa.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/msaa.dir/msaa.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests/msaa.c -o CMakeFiles/msaa.dir/msaa.c.s

tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: tests/CMakeFiles/msaa.dir/flags.make
tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: tests/CMakeFiles/msaa.dir/includes_C.rsp
tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: ../deps/getopt.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/msaa.dir/__/deps/getopt.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c

tests/CMakeFiles/msaa.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/msaa.dir/__/deps/getopt.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c > CMakeFiles/msaa.dir/__/deps/getopt.c.i

tests/CMakeFiles/msaa.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/msaa.dir/__/deps/getopt.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/getopt.c -o CMakeFiles/msaa.dir/__/deps/getopt.c.s

tests/CMakeFiles/msaa.dir/__/deps/glad.c.obj: tests/CMakeFiles/msaa.dir/flags.make
tests/CMakeFiles/msaa.dir/__/deps/glad.c.obj: tests/CMakeFiles/msaa.dir/includes_C.rsp
tests/CMakeFiles/msaa.dir/__/deps/glad.c.obj: ../deps/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object tests/CMakeFiles/msaa.dir/__/deps/glad.c.obj"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/msaa.dir/__/deps/glad.c.obj   -c /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c

tests/CMakeFiles/msaa.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/msaa.dir/__/deps/glad.c.i"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c > CMakeFiles/msaa.dir/__/deps/glad.c.i

tests/CMakeFiles/msaa.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/msaa.dir/__/deps/glad.c.s"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/deps/glad.c -o CMakeFiles/msaa.dir/__/deps/glad.c.s

# Object files for target msaa
msaa_OBJECTS = \
"CMakeFiles/msaa.dir/msaa.c.obj" \
"CMakeFiles/msaa.dir/__/deps/getopt.c.obj" \
"CMakeFiles/msaa.dir/__/deps/glad.c.obj"

# External object files for target msaa
msaa_EXTERNAL_OBJECTS =

tests/msaa.exe: tests/CMakeFiles/msaa.dir/msaa.c.obj
tests/msaa.exe: tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj
tests/msaa.exe: tests/CMakeFiles/msaa.dir/__/deps/glad.c.obj
tests/msaa.exe: tests/CMakeFiles/msaa.dir/build.make
tests/msaa.exe: src/libglfw3dll.a
tests/msaa.exe: tests/CMakeFiles/msaa.dir/linklibs.rsp
tests/msaa.exe: tests/CMakeFiles/msaa.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable msaa.exe"
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && "/C/Program Files/CMake/bin/cmake.exe" -E remove -f CMakeFiles/msaa.dir/objects.a
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/ar.exe cr CMakeFiles/msaa.dir/objects.a @CMakeFiles/msaa.dir/objects1.rsp
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && /C/MinGW/bin/gcc.exe    -Wl,--whole-archive CMakeFiles/msaa.dir/objects.a -Wl,--no-whole-archive  -o msaa.exe -Wl,--out-implib,libmsaa.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/msaa.dir/linklibs.rsp

# Rule to build all files generated by this target.
tests/CMakeFiles/msaa.dir/build: tests/msaa.exe

.PHONY : tests/CMakeFiles/msaa.dir/build

tests/CMakeFiles/msaa.dir/clean:
	cd /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/msaa.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/msaa.dir/clean

tests/CMakeFiles/msaa.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MSYS Makefiles" /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1 /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests /C/MinGW/msys/1.0/home/Makayla/glfw-3.2.1/build/tests/CMakeFiles/msaa.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/msaa.dir/depend

