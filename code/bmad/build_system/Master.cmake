#-----------------------------------------------------------
# Master CMake file.
# Implements the ACC build system.
# How to include:
#      include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
# Found in CMakeLists.txt files in project directories.

#-----------------------------------------------------------
# Set link_directories relative path composition to use new
# behvaior that appends relative path information to the
# CMAKE_CURRENT_SOURCE_DIR value.
#-----------------------------------------------------------
cmake_policy (SET CMP0015 NEW)


#-----------------------------------------------------------
# Disable the inclusion of RPATH from the Acc Build System
#-----------------------------------------------------------
SET (CMAKE_SKIP_RPATH TRUE)

#-----------------------------------------------------------
# Set CESR_FLAGS depening on OS type
#-----------------------------------------------------------
IF (${WIN32})
    SET (CESR_FLAGS "-DCESR_WINCVF")
ELSE ()
    SET (CESR_FLAGS "-DCESR_UNIX -DCESR_LINUX")
ENDIF ()

#------------------------------------------
# Honor requests for compiling with openmp
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_OPENMP})
  SET (ACC_ENABLE_OPENMP 1)
ELSE ()
  SET (ACC_ENABLE_OPENMP 0)
ENDIF ()


#------------------------------------------
# Honor requests for compiling with OpenMPI
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_MPI})
  SET (ACC_ENABLE_MPI 1)
ELSE ()
  SET (ACC_ENABLE_MPI 0)
ENDIF ()


#------------------------------------------
# Honor requests for gfortran compiling with 
# -O2 flag made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_GFORTRAN_OPTIMIZATION}) 
  IF ("$ENV{DIST_F90}" MATCHES "gfortran")
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG "-O2")
  ELSEIF ("$ENV{ACC_SET_F_COMPILER}" MATCHES "gfortran")
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG "-O2")
  ELSE ()
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG)
  ENDIF ()
ENDIF ()


#-------------------------------------------------------
# Import environment variables that influence the build
#-------------------------------------------------------
set (DISTRIBUTION_BUILD $ENV{DIST_BUILD})

IF (${DISTRIBUTION_BUILD})
  IF ("$ENV{ACC_ENABLE_SHARED}" MATCHES "Y")
    SET (CMAKE_SKIP_INSTALL_RPATH TRUE)
  ENDIF ()
  set (FORTRAN_COMPILER $ENV{DIST_F90})
  set (RELEASE_DIR $ENV{DIST_BASE_DIR})
  set (PACKAGES_DIR ${RELEASE_DIR})

  # Explicitly remove 32-bit Libraries from Linux build PATH for 64-bit builds - RT#43178
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND NOT "$ENV{ACC_FORCE_32_BIT}" MATCHES "Y")
    SET (CMAKE_IGNORE_PATH /usr/lib)
    SET (CMAKE_SYSTEM_IGNORE_PATH /usr/lib)
  ENDIF ()

ELSEIF ("$ENV{ACC_SET_F_COMPILER}" MATCHES "gfortran")
  set (FORTRAN_COMPILER "gfortran")
  set (RELEASE_DIR $ENV{ACC_RELEASE_DIR})
  set (PACKAGES_DIR ${RELEASE_DIR}/packages)
ELSE ()
  set (FORTRAN_COMPILER "ifort")
  set (RELEASE_DIR $ENV{ACC_RELEASE_DIR})
  set (PACKAGES_DIR ${RELEASE_DIR}/packages)
ENDIF ()
  
IF (FORTRAN_COMPILER MATCHES "gfortran")

  SET (RELEASE_NAME $ENV{DIST_BASE_DIR})
  IF ("$ENV{DIST_F90}" MATCHES "gfortran")
    SET (RELEASE_NAME_TRUE "Off-site Distribution")
  ENDIF ()
  SET (COMPILER_CHOICE ${FORTRAN_COMPILER})
  SET (CMAKE_Fortran_COMPILER gfortran)
     IF ("${ACC_ENABLE_OPENMP}")
       SET (COMPILER_SPECIFIC_F_FLAGS "-cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none -fopenmp")
       SET (OPENMP_LINK_LIBS "gomp")
     ELSE ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none")
     ENDIF ()
  SET (COMPILER_SPECIFIC_DEBUG_F_FLAGS "-O0 -fno-range-check -fbounds-check -Wuninitialized")

ELSE ()

  SET (RELEASE_NAME $ENV{ACC_RELEASE})
  SET (RELEASE_NAME_TRUE $ENV{ACC_TRUE_RELEASE})
  SET (CMAKE_Fortran_COMPILER ifort)
     IF ("${ACC_ENABLE_OPENMP}")
       SET (OPENMP_LINK_LIBS "-liomp5")
       EXEC_PROGRAM (ifort ARGS -v OUTPUT_VARIABLE INTEL_VERSION_OUPUT)
       IF (( ${INTEL_VERSION_OUPUT} MATCHES "16" ) OR ( ${INTEL_VERSION_OUPUT} MATCHES "17" ))
	 SET (IFORT_OPENMP_FLAG "-qopenmp")
       ELSE()
	 SET (IFORT_OPENMP_FLAG "-openmp")
       ENDIF ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-fpp ${IFORT_OPENMP_FLAG}")
#       SET (ACC_LINK_FLAGS "${IFORT_OPENMP_FLAG}")
     ELSE ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-fpp")
     ENDIF ()
  SET (COMPILER_SPECIFIC_DEBUG_F_FLAGS "-check bounds -check format -check uninit -warn declarations -ftrapuv")

ENDIF ()

IF (${ACC_ENABLE_MPI})
  SET (CMAKE_Fortran_COMPILER mpifort)
  EXEC_PROGRAM (mpifort ARGS --showme:compile OUTPUT_VARIABLE MPI_COMPILE_FLAGS)
  EXEC_PROGRAM (mpifort ARGS --showme:link OUTPUT_VARIABLE MPI_LINK_FLAGS)
  EXEC_PROGRAM (mpifort ARGS --showme:incdirs OUTPUT_VARIABLE MPI_INC_DIR)
  EXEC_PROGRAM (mpifort ARGS --showme:libdirs OUTPUT_VARIABLE MPI_LIB_DIR)
  EXEC_PROGRAM (mpifort ARGS --showme:libs OUTPUT_VARIABLE MPI_LIBS)
ENDIF ()

#----------------------------------------------------------------
# If any pre-build script is specified, run it before building
# any code.  The pre-build script may generate code or header
# files.
#----------------------------------------------------------------
IF (PREBUILD_ACTION)
  message("Executing pre-build action...")
  EXECUTE_PROCESS (COMMAND ${PREBUILD_ACTION}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
ENDIF ()


get_filename_component(SHORT_DIRNAME ${CMAKE_SOURCE_DIR} NAME)
#-------------------------------------------------------
# Check to see if shared object creation is enabled.
# If enabled, shared object libraries will only be
# created for projects that define 
#   CREATE_SHARED
# in their CmakeLists.txt file.
#
# TWO variables need to be set to produce a shared
# object library for a given project.
#  Shell: ACC_ENABLE_SHARED
#  Cmake: CREATE_SHARED
#-------------------------------------------------------
set(ENABLE_SHARED $ENV{ACC_ENABLE_SHARED})

IF (${LIBNAME})
  project(${LIBNAME})
ENDIF ()


#-----------------------------------
#   System library locators
# (Also provides include directory
# locations.)
#-----------------------------------
find_package(X11)


#-----------------------------------
# C / C++ Compiler flags
#-----------------------------------

IF (CMAKE_SYSTEM_NAME MATCHES "HARDWARE-DEVEL")
  SET (BASE_C_FLAGS)
  SET (BASE_CXX_FLAGS)
ELSE ()
  SET (BASE_C_FLAGS "-Df2cFortran -O0 -std=gnu99 ${CESR_FLAGS} -D_POSIX -D_REENTRANT -Wall -fPIC -Wno-trigraphs -Wno-unused ${ACC_MPI_COMPILER_FLAGS}")
  SET (BASE_CXX_FLAGS "-O0 -Wno-deprecated ${CESR_FLAGS} -D_POSIX -D_REENTRANT -Wall -fPIC -Wno-trigraphs -Wno-unused ${ACC_MPI_COMPILER_FLAGS}")
ENDIF ()

#-----------------------------------                                                                                
# For non-Linux or non-ifort or 
# non-64-bit builds, do not use
# unspported flag option 
# "-mcmodel=medium"                                                                      
#-----------------------------------                                                                                
IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${FORTRAN_COMPILER} MATCHES "ifort" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET (BASE_C_FLAGS "${BASE_C_FLAGS} -mcmodel=medium")
  SET (BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -mcmodel=medium")
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "HARDWARE-DEVEL")
  SET (BASE_C_FLAGS)
  SET (BASE_CXX_FLAGS)
ENDIF ()

#-----------------------------------
# Readline Library Definitions 
#-----------------------------------
SET (READLINE_LINK_LIBS readline)
SET (READLINE_LINK_FLAGS "-lreadline")

#-----------------------------------
# Plotting Library Linker flags
#-----------------------------------

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

IF ($ENV{ACC_PLOT_PACKAGE} MATCHES "plplot")
  SET (PLOT_LIBRARY_F_FLAG "-DCESR_PLPLOT")

  IF (${DISTRIBUTION_BUILD})
    IF (NOT "$ENV{ACC_PLOT_DISPLAY_TYPE}" MATCHES "QT")
      SET (PLOT_LINK_FLAGS "-lX11 $ENV{PLOT_LINK_FLAGS}")
    ELSE ()
      SET (PLOT_LINK_FLAGS "$ENV{PLOT_LINK_FLAGS}")
    ENDIF ()
  ELSE ()
    SET (PLOT_LINK_FLAGS "-lX11 $ENV{PLOT_LINK_FLAGS}")
  ENDIF ()

  IF (${MSYS})
    # Assuming Qt backend:
    SET (SHARED_LINK_LIBS QtSvg4 QtGui4 QtCore4 ${SHARED_LINK_LIBS})
  ELSE ()
    SET (SHARED_LINK_LIBS cairo pango-1.0 pangocairo-1.0 gobject-2.0 ${SHARED_LINK_LIBS})
  ENDIF ()

  SET (ACC_PLOT_INC_DIRS)

ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "pgplot")
  SET (PLOT_LIBRARY_F_FLAG "-DCESR_PGPLOT")
  SET (PLOT_LINK_FLAGS "-lX11 $ENV{PLOT_LINK_FLAGS}")

ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "none")
  SET (PLOT_LIBRARY_F_FLAG "-DCESR_NOPLOT")
  SET (PLOT_LINK_LIBS)
  SET (PLOT_LINK_FLAGS)

ENDIF ()

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND NOT $ENV{ACC_PLOT_PACKAGE} MATCHES "none")
  SET (ACC_PLOT_LIB_DIRS /usr/lib64)
ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET (ACC_PLOT_LIB_DIRS /opt/local/lib /opt/X11/lib)
ENDIF ()

IF (${CMAKE_Fortran_COMPILER} MATCHES "ifort" AND "$ENV{ACC_ENABLE_SHARED}" MATCHES "Y")
  SET (PLOT_LINK_FLAGS "${PLOT_LINK_FLAGS} -lifcore -lifport -limf -lsvml -lintlc")
ENDIF ()


#--------------------------------------
# Fortran Compiler flags
#--------------------------------------
enable_language( Fortran )

SET (BASE_Fortran_FLAGS "-Df2cFortran ${CESR_FLAGS} -u -traceback ${COMPILER_SPECIFIC_F_FLAGS} ${PLOT_LIBRARY_F_FLAG} ${MPI_COMPILE_FLAGS}")

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${FORTRAN_COMPILER} MATCHES "ifort" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -mcmodel=medium")
ENDIF ()

IF ($ENV{ACC_ENABLE_FPIC})
   SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -fPIC")
ENDIF ()

IF (${DISTRIBUTION_BUILD})
    SET (ACC_LINK_FLAGS ${ACC_LINK_FLAGS} ${MPI_LINK_FLAGS} ${PLOT_LINK_FLAGS})
ELSE ()
    SET (ACC_LINK_FLAGS "-lreadline -ltermcap -lcurses -lpthread -lstdc++ -lactivemq-cpp" ${ACC_LINK_FLAGS} ${MPI_LINK_FLAGS} ${PLOT_LINK_FLAGS})
ENDIF ()

IF (${MSYS})
    SET (ACC_LINK_FLAGS)
    SET (SHARED_LINK_LIBS ${SHARED_LINK_LIBS} stdc++ readline termcap gdi32 Comdlg32)
ENDIF ()

SET (ACC_INC_DIRS ${ACC_PLOT_INC_DIRS} ${ACC_INC_DIRS} ${MPI_INC_DIRS})
SET (ACC_LIB_DIRS ${ACC_PLOT_LIB_DIRS} ${ACC_LIB_DIRS} ${MPI_LIB_DIRS})


#--------------------------------------
# Honor requests for debug builds 
# composed with any variation in case.
#--------------------------------------
IF (CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
  SET(DEBUG 1)
ELSE ()
  SET(DEBUG 0)
ENDIF ()


#------------------------------------------
# Honor requests for executable building
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_BUILD_EXES})
  SET(BUILD_EXES 1)
ELSE ()
  SET(BUILD_EXES 0)
ENDIF ()


#------------------------------------------
# Honor requests for test executable 
# building made via environment variable.
#------------------------------------------
IF ($ENV{ACC_BUILD_TEST_EXES})
  SET(BUILD_TEST_EXES 1)
ELSE ()
  SET(BUILD_TEST_EXES 0)
ENDIF ()


#-----------------------------------------
# Print some friendly build information
# and according to the build type, set
# the following:
#    C Flags
#    F90 Flags
#-----------------------------------------
message("")
IF (DEBUG)
  message("Build type           : Debug")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../debug)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/debug)
  set (PACKAGES_OUTPUT_BASEDIR ${PACKAGES_DIR}/debug)
  set (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} ${COMPILER_SPECIFIC_DEBUG_F_FLAGS}")
ELSE ()
  message("Build type           : Production")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../production)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/production)
  set (PACKAGES_OUTPUT_BASEDIR ${PACKAGES_DIR}/production)
  set (BASE_C_FLAGS "${BASE_C_FLAGS} -O2")
  set (BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -O2")
  set (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} ${ACC_GFORTRAN_OPTIMIZATION_FLAG}")
ENDIF ()
IF (${CMAKE_SYSTEM_NAME} MATCHES "HARDWARE-DEVEL")
  message("Build Target         : ${CMAKE_SYSTEM_NAME}")
  message("C Compiler           : ${CMAKE_C_COMPILER}")
  message("Fortran Compiler     : ${CMAKE_Fortran_COMPILER}")
  message("Linker               : ${CMAKE_LINKER}")
  set (BASE_C_FLAGS)
ELSE ()
  message("Linking with release : ${RELEASE_NAME} \(${RELEASE_NAME_TRUE}\)")
  message("C Compiler           : ${CMAKE_C_COMPILER}")
  message("Fortran Compiler     : ${CMAKE_Fortran_COMPILER}")
  message("Plotting Libraries   : ${PLOT_LINK_LIBS}")
ENDIF ()
IF (DEFINED SHARED_LINK_LIBS)
message("Shared Libraries     : ${SHARED_LINK_LIBS}")
ENDIF()
IF ($ENV{ACC_ENABLE_OPENMP})
  IF (${FORTRAN_COMPILER} MATCHES "ifort")
      MESSAGE("OpenMP ifort Flag    : ${IFORT_OPENMP_FLAG}")
  ELSE ()
    message("OpenMP gfortran Flag : -fopenmp")
    message("OpenMP Linker Libs   : ${OPENMP_LINK_LIBS}")
  ENDIF()
ELSE ()
  message("OpenMP Support       : Not Enabled")
ENDIF ()
IF ($ENV{ACC_ENABLE_MPI})
  message("MPI Support          : Enabled")
ELSE()
  message("MPI Support          : Not Enabled")
ENDIF()
message("FFLAGS               : ${FFLAGS} ${FCFLAGS}")
message("${FORTRAN_COMPILER} Compiler Flags : ${BASE_Fortran_FLAGS}")
message("${FORTRAN_COMPILER} Linker Flags   : ${ACC_LINK_FLAGS} ${OPENMP_LINK_LIBS}\n")

#-----------------------------------
# Output path definitions
#-----------------------------------
SET (LIBRARY_OUTPUT_PATH ${OUTPUT_BASEDIR}/lib)
SET (EXECUTABLE_OUTPUT_PATH ${OUTPUT_BASEDIR}/bin)
SET (CMAKE_Fortran_MODULE_DIRECTORY ${OUTPUT_BASEDIR}/modules)
SET (INCLUDE_OUTPUT_PATH ${OUTPUT_BASEDIR}/include)

#-------------------------
#   Include directories
#-------------------------
# TODO: Double each include directory entry to search for a local (../<xxx>) path and THEN
#       to a release-based path?
#
# This leaves the possibility that someone may delete the local library, and initiate a build
# that requires that library while leaving the local souce tree and include files intact.
# This new build will then perform the divergent action of linking against the release library
# but extracting constants and other header information from the LOCAL source tree.  

SET (MASTER_INC_DIRS
  ${X11_INCLUDE_DIR}
  ${INC_DIRS}
  ${OUTPUT_BASEDIR}/include
  ${ROOT_INC}
  ${RELEASE_DIR}/include
  ${OUTPUT_BASEDIR}/modules
  ${RELEASE_OUTPUT_BASEDIR}/modules
  ${RELEASE_OUTPUT_BASEDIR}/include
  ${ACC_INC_DIRS}
)

# If not building a distribution, include the include directories which are not part 
# of the distribution.

IF (${DISTRIBUTION_BUILD})
ELSE ()
  SET (MASTER_INC_DIRS
    ${MASTER_INC_DIRS}
    ${PACKAGES_DIR}/${CMAKE_BUILD_TYPE}/include
    ${PACKAGES_DIR}/forest/include
    ${PACKAGES_DIR}/recipes_c-ansi/include
    ${PACKAGES_DIR}/cfortran/include
    ${PACKAGES_DIR}/${CMAKE_BUILD_TYPE}/include/root
    ${PACKAGES_DIR}/${CMAKE_BUILD_TYPE}/include/activemq-cpp-3.7.0
    ${PACKAGES_OUTPUT_BASEDIR}/modules
  )
ENDIF ()

# If building for HARDWARE-DEVEL, remove the include directories which are not needed 

IF (CMAKE_SYSTEM_NAME MATCHES "HARDWARE-DEVEL")
  SET (MASTER_INC_DIRS)
ENDIF ()


#------------------------------------------------------
# Add local include paths to search list if they exist
#------------------------------------------------------
foreach (dir ${INC_DIRS})
  STRING (FIND ${dir} "../" relative)
  STRING (REPLACE "../" "" dirnew ${dir})
  IF (${relative} EQUAL 0)
    LIST (APPEND MASTER_INC_DIRS ${dir})
    LIST (APPEND MASTER_INC_DIRS ${RELEASE_DIR}/${dirnew})
  ELSE ()
    LIST (APPEND MASTER_INC_DIRS ${dir})
  ENDIF ()
endforeach(dir)


LIST (REMOVE_DUPLICATES MASTER_INC_DIRS)
INCLUDE_DIRECTORIES (${MASTER_INC_DIRS})

#--------------------------------------------
# To avoid a Link Lib path "not found" error,
# when a Distribution Build environment is
# has been envoked and in a target project
# directory not within the BMAD_DISTRIBUTION 
# tree, create a empty lib directory.
#--------------------------------------------

IF (NOT EXISTS ${OUTPUT_BASEDIR}/lib)
  FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/lib)
ENDIF ()

#--------------------------------------------
# To avoid an include path "not found" error,
# when a Distribution Build environment is
# has been envoked and in a target project
# directory not within the BMAD_DISTRIBUTION 
# tree, create a empty include directory.
#--------------------------------------------

IF (NOT EXISTS ${OUTPUT_BASEDIR}/include)
  FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/include)
ENDIF ()

#-----------------------------------
# Link directories - order matters
# Lowest level to highest, i.e. in
# order of increasing abstraction.
#-----------------------------------

SET (MASTER_LINK_DIRS
  /usr/lib64
  ${OUTPUT_BASEDIR}/lib
  ${PACKAGES_OUTPUT_BASEDIR}/lib
  ${RELEASE_OUTPUT_BASEDIR}/lib
  ${ACC_LIB_DIRS}
)

IF (DEBUG)
  IF (IS_DIRECTORY "${PACKAGES_DIR}/debug/lib/root")
    LIST (APPEND MASTER_LINK_DIRS "${PACKAGES_DIR}/debug/lib/root")
  ENDIF ()
ELSE ()
  IF (IS_DIRECTORY "${PACKAGES_DIR}/production/lib/root")
    LIST (APPEND MASTER_LINK_DIRS "${PACKAGES_DIR}/production/lib/root")
  ENDIF ()
ENDIF ()

IF (CMAKE_SYSTEM_NAME MATCHES "HARDWARE-DEVEL")
  SET (MASTER_LINK_DIRS)
ENDIF ()

LINK_DIRECTORIES (${MASTER_LINK_DIRS})


#-------------------------------------------
# Collect list of all source files for all
# supported languages from all directories
# mentioned in project CMakeLists.txt file.
#-------------------------------------------
foreach(dir ${SRC_DIRS})

    set(temp_list)
    file(GLOB temp_list ${dir}/*.c)
    LIST(APPEND c_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.C)
    LIST(APPEND c_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for 
    # all C source files.
    #-----------------------------------
    foreach (file ${c_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
    endforeach()
    LIST(APPEND sources ${c_sources})


    set(temp_list)
    file(GLOB temp_list ${dir}/*.cpp)
    LIST(APPEND cpp_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.cc)
    LIST(APPEND cpp_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.cxx)
    LIST(APPEND cpp_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for 
    # all C source files.
    #-----------------------------------
    foreach (file ${cpp_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CPPFLAGS}")
    endforeach()
    LIST(APPEND sources ${cpp_sources})


    set(temp_list)
    file(GLOB temp_list ${dir}/*.f)
    LIST(APPEND f_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.F)
    LIST(APPEND f_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.f90)
    LIST(APPEND f_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for
    # all Fortran source files.
    #-----------------------------------
    foreach (file ${f_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
    endforeach()
    LIST(APPEND sources ${f_sources})

endforeach(dir)


set(DEPS)


set(TARGETS)

IF (LIBNAME)
  add_library( ${LIBNAME} STATIC ${sources} )
  LIST(APPEND TARGETS ${LIBNAME})
  SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES OUTPUT_NAME ${LIBNAME})
ENDIF ()


#----------------------------------------------------------------
# Copy the contents of a config directory to 
# <DIR>/../config/${LIBNAME} if one exists.
#----------------------------------------------------------------
IF (IS_DIRECTORY "../config")
  message("Copying config directory contents to ${CMAKE_SOURCE_DIR}/../config/${SHORT_DIRNAME}...")
  file (MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/../config")
  EXECUTE_PROCESS (COMMAND cp -rfu ${CMAKE_SOURCE_DIR}/config/. ${CMAKE_SOURCE_DIR}/../config/${SHORT_DIRNAME})
ENDIF ()


#----------------------------------------------------------------
# For selectively producing shared object libraries (.so files).
#
# Shell variable ACC_ENABLE_SHARED needs to be set to 
#  "Y" or "true" or "1"
#    - AND - 
# set (CREATE_SHARED true) needs to be present in the individual
#  project's CMakeLists.txt file.
# 
# Now works correctly with gmake -j values greater than 1
#----------------------------------------------------------------
IF (DEFINED SHARED_LINK_LIBS)
MESSAGE ("SHARED DEPS          : ${SHARED_DEPS}\n")
ENDIF ()
IF (ENABLE_SHARED AND CREATE_SHARED)
  ADD_LIBRARY (${LIBNAME}-shared SHARED ${sources})
  ADD_DEPENDENCIES (${LIBNAME}-shared ${LIBNAME}) 
  LIST (APPEND TARGETS ${LIBNAME}-shared)
  SET_TARGET_PROPERTIES (${LIBNAME}-shared PROPERTIES OUTPUT_NAME ${LIBNAME})
  TARGET_LINK_LIBRARIES (${LIBNAME}-shared ${SHARED_DEPS} ${SHARED_LINK_LIBS})
ENDIF ()


#---------------------------------------------------------------
# Add each TEST EXE build description file mentioned in the
# project's CMakeLists.txt file, as requested.
#---------------------------------------------------------------
IF (BUILD_TEST_EXES)
  SET (EXE_SPECS ${TEST_EXE_SPECS} ${EXE_SPECS})
ENDIF ()


#---------------------------------------------------------------
# Process each EXE build description file mentioned in the
# project's CMakeLists.txt file.
#---------------------------------------------------------------
foreach(exespec ${EXE_SPECS})

  set(CFLAGS)
  set(FFLAGS)
  set(CPPFLAGS)

  set(SRC_DIRS)
  set(c_sources)
  set(cpp_sources)
  set(f_sources)

  include(${exespec})


  # TODO: Convert this to macro or function?
  foreach(dir ${INC_DIRS})
    STRING(FIND ${dir} "../" relative)
    STRING(REPLACE "../" "" dirnew ${dir})
    IF (${relative} EQUAL 0)
      LIST(APPEND MASTER_INC_DIRS ${dir})
      LIST(APPEND MASTER_INC_DIRS ${RELEASE_DIR}/${dirnew})
    ELSE ()
      LIST(APPEND MASTER_INC_DIRS ${dir})
    ENDIF ()
  endforeach(dir)
  LIST(REMOVE_DUPLICATES MASTER_INC_DIRS)
  INCLUDE_DIRECTORIES(${MASTER_INC_DIRS})

  set(DEPS ${LINK_LIBS})

  #----------------------------------------------------------------
  # Make this project's EXE build depend upon the product of each
  # build that is listed as a dependency.  If those binaries
  # are newer than this project's EXE, relink this project's EXE.
  # Only invoke add_library to tie in external dependencies a
  # single time for each unique target.
  #----------------------------------------------------------------
  foreach(dep ${DEPS})

    IF(${LIBNAME} MATCHES ${dep})
    ELSE()
      LIST(FIND TARGETS ${dep} DEP_SEEN)
      IF(${DEP_SEEN} EQUAL -1)
        IF (EXISTS ${OUTPUT_BASEDIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${OUTPUT_BASEDIR}/lib/lib${dep}.a)
        ELSEIF (EXISTS ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.a)
        ELSEIF (EXISTS ${RELEASE_DIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${RELEASE_DIR}/lib/lib${dep}.a)
        ENDIF ()
      ENDIF()
    ENDIF()

  endforeach(dep)

  LINK_DIRECTORIES( ${LINK_DIRS} )

  #----------------------------------------------
  # Honor propagated control variable to build
  # any EXEs provided in a list from the main
  # CMakeLists.txt file.
  #----------------------------------------------
  IF (BUILD_EXES)
    set(BUILD_EXE_TOGGLE "")
  ELSE()
    set(BUILD_EXE_TOGGLE "EXCLUDE_FROM_ALL")
  ENDIF ()


  #-------------------------------------------
  # Collect list of all source files for all
  # supported languages from all directories
  # mentioned in project CMakeLists.txt file.
  #-------------------------------------------
  foreach(dir ${SRC_DIRS})

      set(temp_list)
      file(GLOB temp_list ${dir}/*.c)
      LIST(APPEND c_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.C)
      LIST(APPEND c_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for 
      # all C source files.
      #-----------------------------------
      foreach (file ${c_sources})
        set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${c_sources})


      set(temp_list)
      file(GLOB temp_list ${dir}/*.cpp)
      LIST(APPEND cpp_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.cc)
      LIST(APPEND cpp_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.cxx)
      LIST(APPEND cpp_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for 
      # all C source files.
      #-----------------------------------
      foreach (file ${cpp_sources})
	set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CPPFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${cpp_sources})


      set(temp_list)
      file(GLOB temp_list ${dir}/*.f)
      LIST(APPEND f_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.F)
      LIST(APPEND f_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.f90)
      LIST(APPEND f_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for
      # all Fortran source files.
      #-----------------------------------
      foreach (file ${f_sources})
	set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${f_sources})

  endforeach(dir)


  #----------------------------------------------
  # Apply user-specified compiler flags to each 
  # file being built into the executable.
  #----------------------------------------------
  SET(COMPILER_FLAGS "${BASE_C_FLAGS} ${COMPILER_FLAGS} ${CFLAGS}")
  SET(COMPILER_FLAGS "${BASE_Fortran_FLAGS} ${COMPILER_FLAGS} ${FFLAGS}")
  if (COMPILER_FLAGS)
    foreach(srcfile ${SRC_FILES})
      set(ext_match)
      STRING(FIND ${srcfile} ".c" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".C" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cpp" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cc" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cxx" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
      ENDIF()
      #------------
      set(ext_match)
      STRING(FIND ${srcfile} ".f90" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".f" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".F" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
    endforeach(srcfile)
  endif ()


  # Actually request creation of executable target
  add_executable(${EXENAME}-exe
      ${BUILD_EXE_TOGGLE}
      ${SRC_FILES}
  )

  SET_TARGET_PROPERTIES(${EXENAME}-exe
          PROPERTIES
          OUTPUT_NAME
          ${EXENAME}
  )

  IF (DEFINED LINKER_LANGUAGE_PROP)
    SET_TARGET_PROPERTIES (${EXENAME}-exe PROPERTIES LINKER_LANGUAGE ${LINKER_LANGUAGE_PROP})
    UNSET (LINKER_LANGUAGE_PROP)
  ENDIF ()

  #----------------------------------------------
  # When linking an executable using differnet 
  # linker then the one used to build libraries 
  # in the LINK_LIBS statement, the user must 
  # set IMPLICIT_LINK_LIBS to the <lang> of the  
  # libraries in LINK_LIBS. RT#37678 
  #----------------------------------------------
  IF (DEFINED IMPLICIT_LINK_LIBS)
    SET (IMPLICIT_LINKER_LIBRARIES ${CMAKE_${IMPLICIT_LINK_LIBS}_IMPLICIT_LINK_LIBRARIES})
    MESSAGE ("Implicit ${FORTRAN_COMPILER} Linker Flags   : ${IMPLICIT_LINKER_LIBRARIES}\n")
    UNSET (IMPLICIT_LINK_LIBS)
  ENDIF ()

  # Create map file output directory if it doesn't yet exist.
  IF (IS_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ELSE ()
    FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ENDIF ()

  # Set up linking for the executable.
  # Always produce a map file.  It is placed in the ../<build_type>/map directory
  # created during build setup.
  set (MAPLINE "-Wl,-Map=${OUTPUT_BASEDIR}/map/${EXENAME}.map")
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET (MAPLINE "-Wl,-map -Wl,${OUTPUT_BASEDIR}/map/${EXENAME}.map")
  ENDIF ()

  set(STATIC_FLAG "")
  set(SHARED_FLAG "")
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    IF (FORTRAN_COMPILER MATCHES "ifort")
      set(STATIC_FLAG "-Wl,-Bstatic")
      set(SHARED_FLAG "-Wl,-Bdynamic")
    ENDIF ()
  ELSEIF (${MSYS})
      # Link all libraries statically:
      SET (STATIC_FLAG "-static")
      SET (SHARED_FLAG "")
  ENDIF ()

  IF (ENABLE_SHARED AND CREATE_SHARED)
  ELSE()
    set (EXTRA_SHARED_LINK_LIBS "")
  ENDIF ()

  IF (CMAKE_SYSTEM_NAME MATCHES "HARDWARE-DEVEL")
    foreach(file ${SRC_FILES})
      set(OBJ_FILES "${OBJ_FILES} CMakeFiles/${EXENAME}-exe.dir/${file}.o")
    endforeach(file)
    set (CMAKE_C_LINK_EXECUTABLE "${CMAKE_LINKER} ${LINK_FLAGS} ${OBJ_FILES} -o ${OUTPUT_BASEDIR}/bin/${EXENAME}.exe")
    TARGET_LINK_LIBRARIES(${EXENAME}-exe
      )
  ELSE ()
    TARGET_LINK_LIBRARIES(${EXENAME}-exe
      ${STATIC_FLAG} ${LINK_LIBS} 
      ${SHARED_FLAG} ${SHARED_LINK_LIBS} ${EXTRA_SHARED_LINK_LIBS}
      ${X11_LIBRARIES} ${ACC_LINK_FLAGS} ${OPENMP_LINK_LIBS}
      ${LINK_FLAGS} ${MAPLINE} ${IMPLICIT_LINKER_LIBRARIES}
      )
  ENDIF ()

  SET(CFLAGS)
  SET(FFLAGS)
  SET(COMPILER_FLAGS)
  SET(LINK_FLAGS)

  # Copy all header files from the project's specified include directory into the output include directory.
  SET (INCLUDE_INSTALL_DIR ${INCLUDE_OUTPUT_PATH})
  FOREACH (inc_dir ${DIR_OF_INCLUDES_TO_MOVE})
    FILE (GLOB inc_files "${inc_dir}/*.h")
    FILE (COPY ${inc_files} DESTINATION ${INCLUDE_INSTALL_DIR}) 
  ENDFOREACH (inc_dir)

endforeach(exespec)


#-------------------------------------------------------------------
# If a shared object library has been built, and if a Makefile.mex
# file exists for building Matlab MEX wrappers, call that makefile.
#-------------------------------------------------------------------
FOREACH(target ${TARGETS})

  IF(target MATCHES ${LIBNAME}-shared)
    IF(EXISTS ../Makefile.mex)
      ADD_CUSTOM_COMMAND(TARGET ${LIBNAME}-shared
      POST_BUILD
      COMMAND gmake -f Makefile.mex
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      )
    ENDIF(EXISTS ../Makefile.mex)
  ENDIF()

ENDFOREACH()
