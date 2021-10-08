// (C) Copyright Daan Camps, Mercy Amankwah, E. Wes Bethel, Talita Perciano and
//               Roel Van Beeumen

#ifndef qpixl_frqi_util_hpp
#define qpixl_frqi_util_hpp

#include "qclab/util.hpp"
#include "qpixl/util.hpp"
#include <cmath>
#include <array>
#include <complex>

/// QPIXL::FRQI namespace
namespace qpixl::frqi {


  /**
   * \brief Computes a scaled fast Walsh-Hadamard transform in binary ordering.
   *
   * The size of the input vector should be a power of 2.
   */
  template <typename V>
  void sfwht( V& a ) {
    using T = typename V::value_type ;
    size_t n = a.size() ;
    assert( qpixl::ispow2( n ) ) ;
    size_t k = qpixl::ilog2( n ) ;

    for ( size_t j = 1; j < n; j*=2 ) {
      #pragma omp parallel for
      for ( size_t i = 0; i < n; i++ ) {
        if ( ( i & j ) == 0 ) {
          const size_t j1 = i + j ;
          const T x = a[i] ;
          const T y = a[j1] ;

          a[i]  = ( x + y ) / T(2) ;
          a[j1] = ( x - y ) / T(2) ;
        }
      }
    }
  }
  
  /**
   * \brief In-place computation of an inverse scaled fast Walsh-Hadamard 
   * transform in binary order of an input vector `a`.
   * 
   * The size of the input vector should be a power of 2.
   */
  template <typename V>
  void isfwht( V& a ) {
    using T = typename V::value_type ;
    size_t n = a.size() ;
    assert( qpixl::ispow2( n ) ) ;
    size_t k = qpixl::ilog2( n ) ;

    for ( size_t j = 1; j < n; j*=2 ) {
      #pragma omp parallel for
      for ( size_t i = 0; i < n; i++ ) {
        if ( ( i & j ) == 0 ) {
          const size_t j1 = i + j ;
          const T x = a[i] ;

          a[i]  = ( a[i] + a[j1] ) ;
          a[j1] = ( x    - a[j1] ) ;
        }
      }
    }
  }

  /// Gray code of x
  static inline size_t grayCode( const size_t x ) {
    return x ^ ( x>>1 );
  }


  /// \brief Permutes `a` by Gray code permutation in `b`
  template <typename V>
  inline void grayPermutation( const V& a, V& b){
    using T = typename V::value_type ;
    size_t n = a.size() ;
    assert( n == b.size() );
    #pragma omp parallel for
    for (size_t k=0; k<n; k++ ){
      b[k] = a[qpixl::frqi::grayCode(k)];
    }
  }

  /// \brief Permutes `a` by an inverse Gray code permutation in `b`
  template <typename V>
  inline void invGrayPermutation( const V& a, V& b){
    using T = typename V::value_type ;
    size_t n = a.size() ;
    assert( n == b.size() );
    #pragma omp parallel for
    for (size_t k=0; k<n; k++ ){
      b[qpixl::frqi::grayCode(k)] = a[k];
    }
  }  

  /// \brief In-place Gray code permutation of `a`
  template <typename V>
  void grayPermutation( V& a ) {
    using T = typename V::value_type ;
    size_t n = a.size() ;
    size_t z = 1, u = 0, v = 0, cl = 1 ;

    for ( size_t ldm = 1, m = 2; m < n; ++ldm, m<<=1 ) {
      z <<= 1 ;
      v <<= 1 ;
      if ( qpixl::ispow2( ldm ) ){
        ++z ;
        cl <<= 1 ;
      }
      else ++v ;

      do{
        u = (u - v) & v ;
        size_t i = z | u ;
        T t = a[i] ;
        size_t g = qpixl::frqi::grayCode(i) ;
        for ( size_t k = cl-1 ; k != 0; --k ){
          a[i] = a[g] ;
          i = g ;
          g = qpixl::frqi::grayCode( i ) ;
        }
        a[i] = t ;
      }
      while ( u ) ;
    }
  }

  /// \brief In-place inverse Gray code permutation of `a`
  template <typename V>
  void invGrayPermutation( V& a ) {
    using T = typename V::value_type ;
    size_t n = a.size() ;
    size_t z = 1, u = 0, v = 0, cl = 1 ;

    for ( size_t ldm = 1, m = 2; m < n; ++ldm, m<<=1 ) {
      z <<= 1 ;
      v <<= 1 ;
      if ( qpixl::ispow2( ldm ) ){
        ++z ;
        cl <<= 1 ;
      }
      else ++v ;

      do{
        u = (u - v) & v ;
        size_t i = z | u ;
        T t = a[i] ;
        size_t g = qpixl::frqi::grayCode(i) ;
        for ( size_t k = cl-1 ; k != 0; --k ){
          T tt = a[g] ;
          a[g] = t ;
          t = tt ;
          g = qpixl::frqi::grayCode( g ) ;
        }
        a[g] = t ;
      }
      while ( u ) ;
    }
  }

  /// \brief convert a vector containing grayscale data to angles for FRQI
  template< typename V>
  void convertToAngles( V& a , const size_t maxval ) {
    using R = typename V::value_type ;
    const R pi2 = 2 * std::atan(1) ;
    const R scal = pi2 / maxval ;
    std::for_each( a.begin(), a.end(), [&]( R& n ){ n *= scal; } ) ;
  }

  /// \brief convert a vector containing angles for FRQI to grayscale data
  template< typename V>
  void convertToGrayscale( V& a , const size_t maxval ) {
    using R = typename V::value_type ;
    const R pi2 = 2 * std::atan(1) ;
    const R scal = maxval / pi2 ;
    std::for_each( a.begin(), a.end(), [&]( R& n ){ n *= scal; } ) ; 
  }
} // namespace qpixl::frqi

#endif