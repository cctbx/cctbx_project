cmake_minimum_required(VERSION 3.10.0)
project(suitename)

if (WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif ()

set(suitename_SOURCES
	suitename.c
	suiteninit.c
	suitenscrt.c
	suiteninpt.c
	suitenout.c
	suitenutil.c
)

set(suitename_HEADERS
	suitename.h
	suiteninit.h
	suitenscrt.h
	suiteninpt.h
	suitenout.h
	suitenutil.h
)

add_executable(suitename ${suitename_SOURCES} ${suitename_HEADERS})
if (NOT WIN32)
  target_link_libraries(suitename m)
endif()
