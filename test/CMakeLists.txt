add_executable( qpixl_tests run.cpp
                            pgm.cpp
                            util.cpp
                            frqi/util.cpp                           
                            frqi/circuit.cpp
              )
target_link_libraries( qpixl_tests PUBLIC qpixlpp qclabpp gtest )
target_include_directories( qpixl_tests PUBLIC ${PROJECT_SOURCE_DIR}/test )

add_executable( qpixl_frqi_time_sfwht_gray frqi/sfwht_gray.cpp )
target_link_libraries( qpixl_frqi_time_sfwht_gray PUBLIC qpixlpp qclabpp )
target_include_directories( qpixl_frqi_time_sfwht_gray PUBLIC ${PROJECT_SOURCE_DIR}/test )