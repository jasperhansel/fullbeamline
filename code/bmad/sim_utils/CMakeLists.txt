cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(LIBNAME sim_utils)

set (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

set(INC_DIRS 
)

set(SRC_DIRS
  file 
  interfaces 
  io 
  math 
  matrix 
  misc
  optimizers
  plot 
  string
  geodesic_lm
)

# set (FFLAGS "-pg")  # For profiling
# set (CFLAGS "-pg")

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET (SHARED_DEPS
    recipes_f-90_LEPP 
    ${PLOT_LINK_LIBS}
    forest 
    fgsl
    gsl
    gslcblas
    lapack95
    lapack
    blas
    X11
    readline
  )
ENDIF ()

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





