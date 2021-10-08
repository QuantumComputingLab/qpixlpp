#include <gtest/gtest.h>
#include "qpixl/pgm.hpp"
#include "qpixl/frqi/util.hpp"

template <typename R>
void test_qpixl_pgm_readPGMA( ) {
  
  std::vector< R > data ;
  std::string fileName = "../test/test.pgm" ;
  size_t nrows, ncols, maxval ;
  int s ;

  // no zero padding
  s = qpixl::readPGMA( fileName , data , nrows , ncols , maxval , 0 ) ;
  EXPECT_EQ( s , 0  ) ;
  EXPECT_EQ( nrows , 600 ) ;
  EXPECT_EQ( ncols , 600 ) ;
  EXPECT_EQ( maxval , 255 ) ;
  EXPECT_EQ( data.size() , 360000 ) ;

  // embed the image
  s = qpixl::readPGMA( fileName , data , nrows , ncols , maxval , 1 ) ;
  EXPECT_EQ( s , 0  ) ;
  EXPECT_EQ( nrows , 600 ) ;
  EXPECT_EQ( ncols , 600 ) ;
  EXPECT_EQ( maxval , 255 ) ;
  EXPECT_EQ( data.size() , 1048576 ) ;

  // embed the vector
  s = qpixl::readPGMA( fileName , data , nrows , ncols , maxval , 2 ) ;
  EXPECT_EQ( s , 0  ) ;
  EXPECT_EQ( nrows , 600 ) ;
  EXPECT_EQ( ncols , 600 ) ;
  EXPECT_EQ( maxval , 255 ) ;
  EXPECT_EQ( data.size() , 524288 ) ;
  
}

template <typename R>
void test_qpixl_pgm() {
  test_qpixl_pgm_readPGMA< R >() ;
}

/*
 * float
 */
TEST( qpixl_pgm , float ) {
  test_qpixl_pgm< float >() ;
}


/*
 * double
 */
TEST( qpixl_pgm , double ) {
  test_qpixl_pgm< double >() ;
}