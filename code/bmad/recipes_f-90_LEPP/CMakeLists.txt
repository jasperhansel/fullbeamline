cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(LIBNAME  recipes_f-90_LEPP )

set(INC_DIRS 
)

set(SRC_DIRS
  lib_src
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
