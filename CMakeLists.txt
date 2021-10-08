cmake_minimum_required( VERSION 3.16 FATAL_ERROR )

# project
project( QPIXLPP VERSION 0.1.0 LANGUAGES CXX )

# fetch content
include( FetchContent )

# qclabpp
FetchContent_Declare( qclabpp
  GIT_REPOSITORY https://github.com/QuantumComputingLab/qclabpp.git
)
FetchContent_MakeAvailable( qclabpp )

# GTest
FetchContent_Declare( gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
)
FetchContent_MakeAvailable( gtest )

# qpixl/include
add_library( qpixlpp INTERFACE )
target_include_directories( qpixlpp INTERFACE ${PROJECT_SOURCE_DIR}/include )
target_compile_features( qpixlpp INTERFACE cxx_std_20 )

# qpixl/test
add_subdirectory( test )

# qpixl/examples
add_subdirectory( examples )

# openmp
find_package( OpenMP )
if( TARGET OpenMP::OpenMP_CXX )
   target_link_libraries( qpixlpp INTERFACE OpenMP::OpenMP_CXX )
endif()

# documentation
find_package( Doxygen )
if( DOXYGEN_FOUND )
  configure_file( ${CMAKE_SOURCE_DIR}/doc/doxygen/doxygen.dox.in
    ${CMAKE_BINARY_DIR}/doxygen.dox @ONLY)
  add_custom_target( qpixl_doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen.dox
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with doxygen" VERBATIM )
endif()

