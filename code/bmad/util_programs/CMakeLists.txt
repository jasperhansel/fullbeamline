set (LIBNAME util_programs)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})


set(SRC_DIRS 
  modules
)

set (EXE_SPECS 
  #  cmake.ansga2
  #  cmake.aspea2
  cmake.envelope_ibs
  cmake.ptc_flat_file_to_bmad
  cmake.make_a_matching_knob
  cmake.wake_fit
  cmake.bmad_to_autocad
  cmake.bmad_to_blender
  cmake.bmad_to_csrtrack
  cmake.bmad_to_gpt
  cmake.bmad_to_mad_sad_and_xsif
  cmake.compare_tracking_methods_plot
  cmake.compare_tracking_methods_text
  cmake.element_attributes
  cmake.f77_to_f90
  cmake.lattice_cleaner
  cmake.sad_to_bmad_postprocess
)

if($ENV{ACC_ENABLE_MPI})
  list(APPEND EXE_SPECS "cmake.moga")
endif()

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
