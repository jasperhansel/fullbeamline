# cmake/modules/drivers-finish.cmake
#
# Finalize driver initializations
#
# Copyright (C) 2006-2015  Alan W. Irwin
#
# This file is part of PLplot.
#
# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published
# by the Free Software Foundation; version 2 of the License.
#
# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with the file PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

# Module to finalize device options.

# Results are contained in the following variables:
# DRIVERS_LIST: list of device drivers (e.g., gd for the png and jpeg devices)
# where at least one of the associated devices is enabled.

# DEVICES_LIST: list of devices (e.g. png, jpeg), where device is enabled.

# FILE_DEVICES_LIST: list of colon-separated information for enabled
# file devices.  The information included is
# <devicename>:<drivername>:<familied>, where familied is ON or OFF
# (depending on whether familying is required for multi-page results).

# INTERACTIVE_DEVICES_LIST: list of colon-separated information for
# enabled interactive devices.  The information included is
# <devicename>:<drivername>.

# N.B. DRIVERS_LIST is used below to create data to help test plug-in
# devices.  Also, it is output in the summary.  DEVICES_LIST is used
# to iterate through all enabled devices (notably in the drivers and
# src subdirectories).  Also, it is output in the summary. Both
# FILE_DEVICES_LIST and INTERACTIVE_DEVICES_LIST are used in the
# examples subdirectory to set up testing of file devices and
# interactive devices.
set(DRIVERS_LIST)
set(DEVICES_LIST)
set(FILE_DEVICES_LIST)
set(INTERACTIVE_DEVICES_LIST)
foreach(DRIVERS_DEVICE ${DRIVERS_DEVICE_LIST})
  string(REGEX REPLACE "^(.*):.*:.*:.*:.*$" "\\1" DEVICE ${DRIVERS_DEVICE})
  string(REGEX REPLACE "^.*:(.*):.*:.*:.*$" "\\1" DRIVER ${DRIVERS_DEVICE})
  if(PLD_${DEVICE})
    list(APPEND DEVICES_LIST ${DEVICE})
    string(REGEX REPLACE "^.*:.*:.*:(.*):.*$" "\\1" KIND ${DRIVERS_DEVICE})
    string(REGEX REPLACE "^.*:.*:.*:.*:(.*)$" "\\1" REQUIRE_FAMILYING ${DRIVERS_DEVICE})
    if(KIND STREQUAL "F")
      list(APPEND FILE_DEVICES_LIST "${DEVICE}:${DRIVER}:${REQUIRE_FAMILYING}")
    elseif(KIND STREQUAL "I")
      list(APPEND INTERACTIVE_DEVICES_LIST ${DEVICE}:${DRIVER})
    endif(KIND STREQUAL "F")

    set(APPEND_DRIVER ON)
    foreach(DRIVER_IN_LIST ${DRIVERS_LIST})
      if(DRIVER STREQUAL "${DRIVER_IN_LIST}")
        set(APPEND_DRIVER OFF)
      endif(DRIVER STREQUAL "${DRIVER_IN_LIST}")
    endforeach(DRIVER_IN_LIST ${DRIVERS_LIST})
    if(APPEND_DRIVER)
      if(ENABLE_DYNDRIVERS)
      	set(${DRIVER}_INFO)
	if(EXISTS ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.driver_info.in)
          file(STRINGS
  	    ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.driver_info.in
	    ${DRIVER}_INFO
	    )
	endif(EXISTS ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.driver_info.in)
      endif(ENABLE_DYNDRIVERS)
      list(APPEND DRIVERS_LIST ${DRIVER})

      # Prepend driver's source code in the drivers directory to
      # anything set specially for the driver in question in the
      # various special CMake modules for drivers included by
      # drivers.cmake before this module (drivers-finish.cmake) is
      # included by drivers.cmake.
      if(DRIVER STREQUAL "wxwidgets")
	if(OLD_WXWIDGETS)
          set(${DRIVER}_SOURCE
            ${CMAKE_SOURCE_DIR}/drivers/deprecated_${DRIVER}.cpp
            ${CMAKE_SOURCE_DIR}/drivers/deprecated_${DRIVER}_app.cpp
            ${CMAKE_SOURCE_DIR}/drivers/deprecated_${DRIVER}_dc.cpp
            ${CMAKE_SOURCE_DIR}/drivers/deprecated_${DRIVER}_gc.cpp
            ${${DRIVER}_SOURCE}
	    )
          if(HAVE_AGG)
            set(${DRIVER}_SOURCE
              ${CMAKE_SOURCE_DIR}/drivers/deprecated_${DRIVER}_agg.cpp
              ${${DRIVER}_SOURCE}
              )
          endif(HAVE_AGG)
	else(OLD_WXWIDGETS)
          set(${DRIVER}_SOURCE
            ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.cpp
	    ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}_comms.cpp
            ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}_dev.cpp
            ${${DRIVER}_SOURCE}
	    )
	endif(OLD_WXWIDGETS)
      elseif(DRIVER STREQUAL "psttf")
        set(${DRIVER}_SOURCE
	        ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.cc
	        ${${DRIVER}_SOURCE}
	        )
	  elseif(DRIVER STREQUAL "qt")
        set(${DRIVER}_SOURCE
	      ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.cpp
	      ${${DRIVER}_SOURCE}
	      )
      else(DRIVER STREQUAL "wxwidgets")
        set(${DRIVER}_SOURCE
	      ${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.c
	      ${${DRIVER}_SOURCE}
	      )
      endif(DRIVER STREQUAL "wxwidgets")
    endif(APPEND_DRIVER)
  endif(PLD_${DEVICE})
endforeach(DRIVERS_DEVICE)

# Calculate driver information and store it in
# ${CMAKE_BINARY_DIR}/drivers/${DRIVER}.driver_info for each driver to be compared
# at run-time with the same information obtained from the actual
# driver plug-in by test-drv-info as a check of the validity of
# that plug-in (and consistency of the driver code with DRIVERS_DEVICE_LIST
# maintained in drivers-init.cmake).
foreach(DRIVERS_DEVICE ${DRIVERS_DEVICE_LIST})
  string(REGEX REPLACE "^(.*):.*:.*:.*:.*$" "\\1" DEVICE ${DRIVERS_DEVICE})
  string(REGEX REPLACE "^.*:(.*):.*:.*:.*$" "\\1" DRIVER ${DRIVERS_DEVICE})
  if(${DRIVER}_INFO)
    if(NOT PLD_${DEVICE})
      set(DEVICE_INFO_MATCHED)
      # Must remove corresponding data from ${DRIVER}_INFO
      #message("DEVICE = ${DEVICE}")
      foreach(DEVICE_INFO ${${DRIVER}_INFO})
        string(REGEX REPLACE "^(.*):.*:.*:.*:.*:.*$" "\\1" DEVICE_INFO_NAME ${DEVICE_INFO})
        #message(STATUS "DEBUG: DEVICE_INFO_NAME = ${DEVICE_INFO_NAME}")
	if(DEVICE STREQUAL "${DEVICE_INFO_NAME}")
	  #There should one and only one match.
	  set(DEVICE_INFO_MATCHED ${DEVICE_INFO})
	endif(DEVICE STREQUAL "${DEVICE_INFO_NAME}")
      endforeach(DEVICE_INFO ${${DRIVER}_INFO})
      if(DEVICE_INFO_MATCHED)
        list(REMOVE_ITEM ${DRIVER}_INFO ${DEVICE_INFO_MATCHED})
      else(DEVICE_INFO_MATCHED)
        message(FATAL_ERROR "${CMAKE_SOURCE_DIR}/drivers/${DRIVER}.driver_info.in not consistent with ${CMAKE_SOURCE_DIR}/cmake/modules/drivers-init.cmake")
      endif(DEVICE_INFO_MATCHED)
    endif(NOT PLD_${DEVICE})
  endif(${DRIVER}_INFO)
endforeach(DRIVERS_DEVICE)

foreach(DRIVERS_DEVICE ${DRIVERS_DEVICE_LIST})
  string(REGEX REPLACE "^.*:(.*):.*:.*:.*$" "\\1" DRIVER ${DRIVERS_DEVICE})
  if(${DRIVER}_INFO)
    file(WRITE ${CMAKE_BINARY_DIR}/drivers/${DRIVER}.driver_info "")
    foreach(DEVICE_INFO ${${DRIVER}_INFO})
      file(APPEND ${CMAKE_BINARY_DIR}/drivers/${DRIVER}.driver_info "${DEVICE_INFO}\n")
    endforeach(DEVICE_INFO ${${DRIVER}_INFO})
  endif(${DRIVER}_INFO)
endforeach(DRIVERS_DEVICE ${DRIVERS_DEVICE_LIST})
