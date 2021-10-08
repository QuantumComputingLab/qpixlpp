#include <gtest/gtest.h>
#include "qpixl/util.hpp"

template <typename T>
void test_qpixl_util_ilog2() {

  EXPECT_EQ( qpixl::ilog2( 1) , 0 ) ;
  EXPECT_EQ( qpixl::ilog2( 2) , 1 ) ;
  EXPECT_EQ( qpixl::ilog2( 4) , 2 ) ;
  EXPECT_EQ( qpixl::ilog2( 8) , 3 ) ;
  EXPECT_EQ( qpixl::ilog2(16) , 4 ) ;
  EXPECT_EQ( qpixl::ilog2(32) , 5 ) ;
  EXPECT_EQ( qpixl::ilog2(64) , 6 ) ;

}

template <typename T>
void test_qpixl_util_ispow2() {
  EXPECT_EQ( qpixl::ispow2(    0), true ) ;
  EXPECT_EQ( qpixl::ispow2(    1), true ) ;
  EXPECT_EQ( qpixl::ispow2(    2), true ) ;
  EXPECT_EQ( qpixl::ispow2(    3), false ) ;
  EXPECT_EQ( qpixl::ispow2( 1024), true ) ;
  EXPECT_EQ( qpixl::ispow2( 1022), false ) ;
}

template <typename T>
void test_qpixl_util_nextpow2() {
  EXPECT_EQ( qpixl::nextpow2(    1),    1 ) ;
  EXPECT_EQ( qpixl::nextpow2(    2),    2 ) ;
  EXPECT_EQ( qpixl::nextpow2(    3),    4 ) ;
  EXPECT_EQ( qpixl::nextpow2(    4),    4 ) ;
  EXPECT_EQ( qpixl::nextpow2(    6),    8 ) ;
  EXPECT_EQ( qpixl::nextpow2(    8),    8 ) ;
  EXPECT_EQ( qpixl::nextpow2(  100),  128 ) ;
  EXPECT_EQ( qpixl::nextpow2( 2000), 2048 ) ;
}

template <typename R>
void test_qpixl_util() {
  test_qpixl_util_ilog2< size_t >() ;
  test_qpixl_util_ispow2< size_t >() ;
  test_qpixl_util_nextpow2< size_t >() ;
}

/*
 * float
 */
TEST( qpixl_util , float ) {
  test_qpixl_util< float >() ;
}


/*
 * double
 */
TEST( qpixl_util , double ) {
  test_qpixl_util< double >() ;
}