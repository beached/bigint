cmake_minimum_required(VERSION 3.13)
project(bigint)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include("${CMAKE_SOURCE_DIR}/glean/CMakeListsCompiler.txt")
include("${CMAKE_SOURCE_DIR}/glean/CMakeLists.txt")

set( HEADER_FOLDER "include" )
set( TEST_FOLDER "tests" )

#find_package(Boost COMPONENTS system iostreams filesystem REQUIRED)
#include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#link_directories(${Boost_LIBRARY_DIRS})

set( HEADER_FILES
	${HEADER_FOLDER}/daw/daw_bigint.h
	${HEADER_FOLDER}/daw/bigint_impl/bigint_impl.h
)

#Tests
enable_testing( )
include_directories("${CMAKE_SOURCE_DIR}/include")
add_custom_target( check )

add_executable(bigint_test EXCLUDE_FROM_ALL tests/bigint.cpp ${HEADER_FILES} )
add_dependencies( check bigint_test )

install( DIRECTORY ${HEADER_FOLDER}/daw DESTINATION include )
