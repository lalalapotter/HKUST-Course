# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build

# Include any dependencies generated for this target.
include ext_build/glfw/src/CMakeFiles/glfw_objects.dir/depend.make

# Include the progress variables for this target.
include ext_build/glfw/src/CMakeFiles/glfw_objects.dir/progress.make

# Include the compile flags for this target's objects.
include ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o: ../ext/glfw/src/context.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/context.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/context.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/context.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/context.c > CMakeFiles/glfw_objects.dir/context.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/context.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/context.c -o CMakeFiles/glfw_objects.dir/context.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o: ../ext/glfw/src/init.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/init.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/init.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/init.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/init.c > CMakeFiles/glfw_objects.dir/init.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/init.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/init.c -o CMakeFiles/glfw_objects.dir/init.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o: ../ext/glfw/src/input.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/input.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/input.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/input.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/input.c > CMakeFiles/glfw_objects.dir/input.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/input.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/input.c -o CMakeFiles/glfw_objects.dir/input.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o: ../ext/glfw/src/monitor.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/monitor.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/monitor.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/monitor.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/monitor.c > CMakeFiles/glfw_objects.dir/monitor.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/monitor.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/monitor.c -o CMakeFiles/glfw_objects.dir/monitor.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.o: ../ext/glfw/src/vulkan.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/vulkan.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/vulkan.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/vulkan.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/vulkan.c > CMakeFiles/glfw_objects.dir/vulkan.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/vulkan.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/vulkan.c -o CMakeFiles/glfw_objects.dir/vulkan.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o: ../ext/glfw/src/window.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/window.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/window.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/window.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/window.c > CMakeFiles/glfw_objects.dir/window.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/window.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/window.c -o CMakeFiles/glfw_objects.dir/window.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.o: ../ext/glfw/src/cocoa_init.m
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/cocoa_init.m.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_init.m

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/cocoa_init.m.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_init.m > CMakeFiles/glfw_objects.dir/cocoa_init.m.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/cocoa_init.m.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_init.m -o CMakeFiles/glfw_objects.dir/cocoa_init.m.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.o: ../ext/glfw/src/cocoa_joystick.m
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/cocoa_joystick.m.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_joystick.m

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/cocoa_joystick.m.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_joystick.m > CMakeFiles/glfw_objects.dir/cocoa_joystick.m.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/cocoa_joystick.m.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_joystick.m -o CMakeFiles/glfw_objects.dir/cocoa_joystick.m.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.o: ../ext/glfw/src/cocoa_monitor.m
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/cocoa_monitor.m.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_monitor.m

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/cocoa_monitor.m.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_monitor.m > CMakeFiles/glfw_objects.dir/cocoa_monitor.m.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/cocoa_monitor.m.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_monitor.m -o CMakeFiles/glfw_objects.dir/cocoa_monitor.m.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.o: ../ext/glfw/src/cocoa_window.m
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/cocoa_window.m.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_window.m

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/cocoa_window.m.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_window.m > CMakeFiles/glfw_objects.dir/cocoa_window.m.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/cocoa_window.m.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_window.m -o CMakeFiles/glfw_objects.dir/cocoa_window.m.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.o: ../ext/glfw/src/cocoa_time.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/cocoa_time.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_time.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/cocoa_time.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_time.c > CMakeFiles/glfw_objects.dir/cocoa_time.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/cocoa_time.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/cocoa_time.c -o CMakeFiles/glfw_objects.dir/cocoa_time.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o: ../ext/glfw/src/posix_tls.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/posix_tls.c.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/posix_tls.c

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/posix_tls.c.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/posix_tls.c > CMakeFiles/glfw_objects.dir/posix_tls.c.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/posix_tls.c.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/posix_tls.c -o CMakeFiles/glfw_objects.dir/posix_tls.c.s

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.o: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/flags.make
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.o: ../ext/glfw/src/nsgl_context.m
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.o"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glfw_objects.dir/nsgl_context.m.o -c /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/nsgl_context.m

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw_objects.dir/nsgl_context.m.i"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/nsgl_context.m > CMakeFiles/glfw_objects.dir/nsgl_context.m.i

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw_objects.dir/nsgl_context.m.s"
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src/nsgl_context.m -o CMakeFiles/glfw_objects.dir/nsgl_context.m.s

glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/context.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/init.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/input.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/monitor.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/vulkan.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/window.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_init.m.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_joystick.m.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_monitor.m.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_window.m.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/cocoa_time.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/posix_tls.c.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/nsgl_context.m.o
glfw_objects: ext_build/glfw/src/CMakeFiles/glfw_objects.dir/build.make

.PHONY : glfw_objects

# Rule to build all files generated by this target.
ext_build/glfw/src/CMakeFiles/glfw_objects.dir/build: glfw_objects

.PHONY : ext_build/glfw/src/CMakeFiles/glfw_objects.dir/build

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/clean:
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw_objects.dir/cmake_clean.cmake
.PHONY : ext_build/glfw/src/CMakeFiles/glfw_objects.dir/clean

ext_build/glfw/src/CMakeFiles/glfw_objects.dir/depend:
	cd /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/ext/glfw/src /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src /Users/zhangcengguang/Documents/GitHub/HKUST-Course/COMP5411/MeshViewer/build/ext_build/glfw/src/CMakeFiles/glfw_objects.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext_build/glfw/src/CMakeFiles/glfw_objects.dir/depend
