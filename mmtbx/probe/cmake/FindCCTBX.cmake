# Copyright (c) 2021, Richardson Lab at Duke
# Licensed under the Apache 2 license

# - Try to find the various packages of CCTBX required to build Probe
# On success this will define the following:
#  CCTBX_FOUND - System has CCTBX
#  CCTBX_INCLUDE_DIRS - The include directories for various CCTBX modules
#  CCTBX_LIBRARIES - The libraries needed to use various CCTBX modules
#  CCTBX_DEFINITIONS - Compiler switches required for using CCTBX

find_package(PkgConfig)

find_path(SCITBX_INCLUDE_DIR scitbx/array_family/operator_traits_builtin.h
	  HINTS
	  	C:/cctbx/build/include
		C:/tmp/cctbx/build/include
		C:/tmp/cctbx_phenix/build/include
	)
find_path(SCITBX_SRC_INCLUDE_DIR scitbx/vec3.h
	  HINTS
	  	C:/cctbx/modules/cctbx_project
		C:/tmp/cctbx/modules/cctbx_project
		C:/tmp/cctbx_phenix/modules/cctbx_project
	)

set(CCTBX_INCLUDE_DIRS ${SCITBX_INCLUDE_DIR} ${SCITBX_SRC_INCLUDE_DIR})

find_library(IOTBX_PDB_LIBRARY
	NAMES
	iotbx_pdb
	PATH_SUFFIXES
	${_libsuffixes}
	HINTS
	"C:/cctbx/build/lib"
	"C:/tmp/cctbx/build/lib"
	"C:/tmp/cctbx_phenix/build/lib"
	PATHS
	C:/usr/local
	/usr/local)

set(CCTBX_LIBRARIES ${IOTBX_PDB_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CCTBX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CCTBX DEFAULT_MSG
                                  IOTBX_PDB_LIBRARY SCITBX_INCLUDE_DIR)

mark_as_advanced(SCITBX_INCLUDE_DIR SCITBX_SRC_INCLUDE_DIR)
