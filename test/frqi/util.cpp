#include <gtest/gtest.h>
#include "qpixl/frqi/util.hpp"

template <typename R>
void test_qpixl_frqi_util_sfwht() {
  const R eps = std::numeric_limits< R >::epsilon() ;

  std::vector< R > a = { 1 , 2 } ;

  {
    qpixl::frqi::sfwht( a );
    EXPECT_EQ( a[0] , R( 1.5) );
    EXPECT_EQ( a[1] , R(-0.5) );
  }

  std::vector< R > b = { 1 , 2 , 3 , 4 };

  {
    qpixl::frqi::sfwht( b );
    EXPECT_NEAR( b[0] , R( 2.5) , eps );
    EXPECT_NEAR( b[1] , R(-0.5) , eps );
    EXPECT_NEAR( b[2] , R(-1.0) , eps );
    EXPECT_NEAR( b[3] , R( 0.0) , eps );
  }
  
  std::vector< R > c = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 };

  {
    qpixl::frqi::sfwht( c );
    EXPECT_NEAR( c[0] , R( 4.5) , eps );
    EXPECT_NEAR( c[1] , R(-0.5) , eps );
    EXPECT_NEAR( c[2] , R(-1.0) , eps );
    EXPECT_NEAR( c[3] , R( 0.0) , eps );
    EXPECT_NEAR( c[4] , R(-2.0) , eps );
    EXPECT_NEAR( c[5] , R( 0.0) , eps );
    EXPECT_NEAR( c[6] , R( 0.0) , eps );
    EXPECT_NEAR( c[7] , R( 0.0) , eps );
  }

  std::vector< R > d = { 1.34 , -3.24 , 2.88 , 16.2 , 0 , -7 , 1.11 , 1.45 ,
                         6.22 , -13.57 , 4.9 , 5.31 , -2.81 , - 3.41 , -1.88 , 
                         10 };
  {
    qpixl::frqi::sfwht( d );
    EXPECT_NEAR( d[ 0] , R(1.09375)  , 10*eps );
    EXPECT_NEAR( d[ 1] , R(0.37625)  , 10*eps );
    EXPECT_NEAR( d[ 2] , R(-3.9025)  , 10*eps );
    EXPECT_NEAR( d[ 3] , R(3.62)     , 10*eps );
    EXPECT_NEAR( d[ 4] , R(1.41125)  , 10*eps );
    EXPECT_NEAR( d[ 5] , R(0.95375)  , 10*eps );
    EXPECT_NEAR( d[ 6] , R(-0.915)   , 10*eps );
    EXPECT_NEAR( d[ 7] , R(1.1425)   , 10*eps );
    EXPECT_NEAR( d[ 8] , R(0.49875)  , 10*eps );
    EXPECT_NEAR( d[ 9] , R(-0.63625) , 10*eps );
    EXPECT_NEAR( d[10] , R(0.085)    , 10*eps );
    EXPECT_NEAR( d[11] , R(-0.465)   , 10*eps );
    EXPECT_NEAR( d[12] , R(1.29125)  , 10*eps );
    EXPECT_NEAR( d[13] , R(-2.87875) , 10*eps );
    EXPECT_NEAR( d[14] , R(-0.5125)  , 10*eps );
    EXPECT_NEAR( d[15] , R(0.1775)   , 10*eps );
  }
}

template <typename R>
void test_qpixl_frqi_util_isfwht() {

  const R eps = std::numeric_limits< R >::epsilon() ;

  std::vector< R > a = { 1.34 , -3.24 , 2.88 , 16.2 , 0 , -7 , 1.11 , 1.45 ,
                           6.22 , -13.57 , 4.9 , 5.31 , -2.81 , - 3.41 , -1.88 ,
                           10 };

  std::vector< R > a_c = { 1.34 , -3.24 , 2.88 , 16.2 , 0 , -7 , 1.11 , 1.45 ,
                           6.22 , -13.57 , 4.9 , 5.31 , -2.81 , - 3.41 , -1.88 ,
                           10 };
  qpixl::frqi::sfwht( a ) ;
  qpixl::frqi::isfwht( a ) ;

  EXPECT_NEAR( a[ 0] , a_c[ 0]  , 100*eps );
  EXPECT_NEAR( a[ 1] , a_c[ 1]  , 100*eps );
  EXPECT_NEAR( a[ 2] , a_c[ 2]  , 100*eps );
  EXPECT_NEAR( a[ 3] , a_c[ 3]  , 100*eps );
  EXPECT_NEAR( a[ 4] , a_c[ 4]  , 100*eps );
  EXPECT_NEAR( a[ 5] , a_c[ 5]  , 100*eps );
  EXPECT_NEAR( a[ 6] , a_c[ 6]  , 100*eps );
  EXPECT_NEAR( a[ 7] , a_c[ 7]  , 100*eps );
  EXPECT_NEAR( a[ 8] , a_c[ 8]  , 100*eps );
  EXPECT_NEAR( a[ 9] , a_c[ 9]  , 100*eps );
  EXPECT_NEAR( a[10] , a_c[10]  , 100*eps );
  EXPECT_NEAR( a[11] , a_c[11]  , 100*eps );
  EXPECT_NEAR( a[12] , a_c[12]  , 100*eps );
  EXPECT_NEAR( a[13] , a_c[13]  , 100*eps );
  EXPECT_NEAR( a[14] , a_c[14]  , 100*eps );
  EXPECT_NEAR( a[15] , a_c[15]  , 100*eps );

}

template <typename R>
void test_qpixl_frqi_util_gray_permutation() {
  std::vector< R > a = { 1 } ;
  std::vector< R > a_out( 1 ) ;

  {
    // with copy
    qpixl::frqi::grayPermutation( a, a_out );
    EXPECT_EQ( a[0] , a_out[0] );
    // in-place
    qpixl::frqi::grayPermutation( a );
    EXPECT_EQ( a[0] , a_out[0] );
  }
  
  std::vector< R > b = { 1 , 2 };
  std::vector< R > b_out( 2 ) ;

  {
    // with copy
    qpixl::frqi::grayPermutation( b, b_out );
    EXPECT_EQ( b[0] , b_out[0] );
    EXPECT_EQ( b[1] , b_out[1] );
    // in-place
    qpixl::frqi::grayPermutation( b );
    EXPECT_EQ( b[0] , b_out[0] );
    EXPECT_EQ( b[1] , b_out[1] );
  } 

  std::vector< R > c = { 1 , 2 , 3 , 4 };
  std::vector< R > c_out( 4 ) ;

  {
    // with copy
    qpixl::frqi::grayPermutation( c, c_out );
    EXPECT_EQ( c[0] , c_out[0] );
    EXPECT_EQ( c[1] , c_out[1] );
    EXPECT_EQ( c[2] , c_out[3] );
    EXPECT_EQ( c[3] , c_out[2] );
    // in-place
    qpixl::frqi::grayPermutation( c );
    EXPECT_EQ( c[0] , c_out[0] );
    EXPECT_EQ( c[1] , c_out[1] );
    EXPECT_EQ( c[2] , c_out[2] );
    EXPECT_EQ( c[3] , c_out[3] );
  } 

  std::vector< R > d = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 };
  std::vector< R > d_out( 8 ) ;

  {
    // with copy
    qpixl::frqi::grayPermutation( d, d_out );
    EXPECT_EQ( d[0] , d_out[0] );
    EXPECT_EQ( d[1] , d_out[1] );
    EXPECT_EQ( d[2] , d_out[3] );
    EXPECT_EQ( d[3] , d_out[2] );
    EXPECT_EQ( d[4] , d_out[7] );
    EXPECT_EQ( d[5] , d_out[6] );
    EXPECT_EQ( d[6] , d_out[4] );
    EXPECT_EQ( d[7] , d_out[5] );
    // in-place
    qpixl::frqi::grayPermutation( d );
    EXPECT_EQ( d[0] , d_out[0] );
    EXPECT_EQ( d[1] , d_out[1] );
    EXPECT_EQ( d[2] , d_out[2] );
    EXPECT_EQ( d[3] , d_out[3] );
    EXPECT_EQ( d[4] , d_out[4] );
    EXPECT_EQ( d[5] , d_out[5] );
    EXPECT_EQ( d[6] , d_out[6] );
    EXPECT_EQ( d[7] , d_out[7] );
  }

  std::vector< R > e = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 
                         10 , 11 , 12 , 13 , 14 , 15 , 16 };
  std::vector< R > e_out( 16 ) ;

  {
    // with copy
    qpixl::frqi::grayPermutation( e, e_out );
    EXPECT_EQ( e[ 0] , e_out[ 0] );
    EXPECT_EQ( e[ 1] , e_out[ 1] );
    EXPECT_EQ( e[ 2] , e_out[ 3] );
    EXPECT_EQ( e[ 3] , e_out[ 2] );
    EXPECT_EQ( e[ 4] , e_out[ 7] );
    EXPECT_EQ( e[ 5] , e_out[ 6] );
    EXPECT_EQ( e[ 6] , e_out[ 4] );
    EXPECT_EQ( e[ 7] , e_out[ 5] );
    EXPECT_EQ( e[ 8] , e_out[15] );
    EXPECT_EQ( e[ 9] , e_out[14] );
    EXPECT_EQ( e[10] , e_out[12] );
    EXPECT_EQ( e[11] , e_out[13] );
    EXPECT_EQ( e[12] , e_out[ 8] );
    EXPECT_EQ( e[13] , e_out[ 9] );
    EXPECT_EQ( e[14] , e_out[11] );
    EXPECT_EQ( e[15] , e_out[10] );
    // in-place
    qpixl::frqi::grayPermutation( e );
    EXPECT_EQ( e[ 0] , e_out[ 0] );
    EXPECT_EQ( e[ 1] , e_out[ 1] );
    EXPECT_EQ( e[ 2] , e_out[ 2] );
    EXPECT_EQ( e[ 3] , e_out[ 3] );
    EXPECT_EQ( e[ 4] , e_out[ 4] );
    EXPECT_EQ( e[ 5] , e_out[ 5] );
    EXPECT_EQ( e[ 6] , e_out[ 6] );
    EXPECT_EQ( e[ 7] , e_out[ 7] );
    EXPECT_EQ( e[ 8] , e_out[ 8] );
    EXPECT_EQ( e[ 9] , e_out[ 9] );
    EXPECT_EQ( e[10] , e_out[10] );
    EXPECT_EQ( e[11] , e_out[11] );
    EXPECT_EQ( e[12] , e_out[12] );
    EXPECT_EQ( e[13] , e_out[13] );
    EXPECT_EQ( e[14] , e_out[14] );
    EXPECT_EQ( e[15] , e_out[15] );
  }
  
}

template <typename R>
void test_qpixl_frqi_util_inv_gray_permutation() {
  std::vector< R > a = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 
                         10 , 11 , 12 , 13 , 14 , 15 , 16  } ;
  std::vector< R > a_c = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 
                         10 , 11 , 12 , 13 , 14 , 15 , 16 };
  std::vector< R > a_out( 16 ) ;
  qpixl::frqi::grayPermutation( a ) ;
    {
    // with copy
    qpixl::frqi::invGrayPermutation( a , a_out ) ;
    EXPECT_EQ( a_out[ 0] , a_c[ 0] );
    EXPECT_EQ( a_out[ 1] , a_c[ 1] );
    EXPECT_EQ( a_out[ 2] , a_c[ 2] );
    EXPECT_EQ( a_out[ 3] , a_c[ 3] );
    EXPECT_EQ( a_out[ 4] , a_c[ 4] );
    EXPECT_EQ( a_out[ 5] , a_c[ 5] );
    EXPECT_EQ( a_out[ 6] , a_c[ 6] );
    EXPECT_EQ( a_out[ 7] , a_c[ 7] );
    EXPECT_EQ( a_out[ 8] , a_c[ 8] );
    EXPECT_EQ( a_out[ 9] , a_c[ 9] );
    EXPECT_EQ( a_out[10] , a_c[10] );
    EXPECT_EQ( a_out[11] , a_c[11] );
    EXPECT_EQ( a_out[12] , a_c[12] );
    EXPECT_EQ( a_out[13] , a_c[13] );
    EXPECT_EQ( a_out[14] , a_c[14] );
    EXPECT_EQ( a_out[15] , a_c[15] );
    // in-place
    qpixl::frqi::invGrayPermutation( a ) ;
    EXPECT_EQ( a[ 0] , a_c[ 0] );
    EXPECT_EQ( a[ 1] , a_c[ 1] );
    EXPECT_EQ( a[ 2] , a_c[ 2] );
    EXPECT_EQ( a[ 3] , a_c[ 3] );
    EXPECT_EQ( a[ 4] , a_c[ 4] );
    EXPECT_EQ( a[ 5] , a_c[ 5] );
    EXPECT_EQ( a[ 6] , a_c[ 6] );
    EXPECT_EQ( a[ 7] , a_c[ 7] );
    EXPECT_EQ( a[ 8] , a_c[ 8] );
    EXPECT_EQ( a[ 9] , a_c[ 9] );
    EXPECT_EQ( a[10] , a_c[10] );
    EXPECT_EQ( a[11] , a_c[11] );
    EXPECT_EQ( a[12] , a_c[12] );
    EXPECT_EQ( a[13] , a_c[13] );
    EXPECT_EQ( a[14] , a_c[14] );
    EXPECT_EQ( a[15] , a_c[15] );
  }

}

template <typename R>
void test_qpixl_frqi_util_convert_angles_grayscale() {
  std::vector< R > a = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 } ;
  std::vector< R > a_c = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 } ;
  size_t maxval = 255 ;

  const R eps = std::numeric_limits< R >::epsilon() ;
  const R pi2 = 2 * std::atan(1) ;
  const R scal = pi2 / maxval ;

  qpixl::frqi::convertToAngles( a , maxval ) ;
  EXPECT_NEAR( a[ 0] , a_c[ 0]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 1] , a_c[ 1]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 2] , a_c[ 2]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 3] , a_c[ 3]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 4] , a_c[ 4]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 5] , a_c[ 5]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 6] , a_c[ 6]*scal  , 10*eps ) ;
  EXPECT_NEAR( a[ 7] , a_c[ 7]*scal  , 10*eps ) ;
  
  qpixl::frqi::convertToGrayscale( a , maxval ) ;
  EXPECT_NEAR( a[ 0] , a_c[ 0]  , 10*eps ) ;
  EXPECT_NEAR( a[ 1] , a_c[ 1]  , 10*eps ) ;
  EXPECT_NEAR( a[ 2] , a_c[ 2]  , 10*eps ) ;
  EXPECT_NEAR( a[ 3] , a_c[ 3]  , 10*eps ) ;
  EXPECT_NEAR( a[ 4] , a_c[ 4]  , 10*eps ) ;
  EXPECT_NEAR( a[ 5] , a_c[ 5]  , 10*eps ) ;
  EXPECT_NEAR( a[ 6] , a_c[ 6]  , 10*eps ) ;
  EXPECT_NEAR( a[ 7] , a_c[ 7]  , 10*eps ) ;

}

template <typename R>
void test_qpixl_frqi_util() {
  test_qpixl_frqi_util_sfwht< R >() ;
  test_qpixl_frqi_util_isfwht< R >() ;
  test_qpixl_frqi_util_gray_permutation< R >() ;
  test_qpixl_frqi_util_inv_gray_permutation< R >() ;
  test_qpixl_frqi_util_convert_angles_grayscale< R >();
}

/*
 * float
 */
TEST( qpixl_frqi_util , float ) {
  test_qpixl_frqi_util< float >() ;
}


/*
 * double
 */
TEST( qpixl_frqi_util , double ) {
  test_qpixl_frqi_util< double >() ;
}

