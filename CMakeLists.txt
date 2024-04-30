# CMake Minimum Required Version
cmake_minimum_required(VERSION 3.22)

# Defining Project
project(BIOS VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED True)

# Defining Variables
set(SOURCES_FOLDER src)
set(LIBRARY_FOLDER libs)
set(
   PATH_INCLUDE
   ${SOURCES_FOLDER}
   ${LIBRARY_FOLDER}
)
set(
   SUB_FOLDERS
   ctrl
   main
   mlib
   pen
   prob 
   sol 
   sur
)

# Adding Executable Target
add_executable(bios ${SOURCES_FOLDER}/main/main.cpp)
target_include_directories(bios PUBLIC ${PATH_INCLUDE})

# Adding Subdirectories
foreach(sub_folder ${SUB_FOLDERS})
   add_subdirectory("${SOURCES_FOLDER}/${sub_folder}")
endforeach()

# Defining Install Configuration
install(
   TARGETS bios
   DESTINATION bin
)