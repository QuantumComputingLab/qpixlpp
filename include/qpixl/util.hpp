// (C) Copyright Daan Camps, Mercy Amankwah, E. Wes Bethel, Talita Perciano and
//               Roel Van Beeumen

#ifndef qpixl_util_hpp
#define qpixl_util_hpp

#include <bit>

/// QPIXL::FRQI namespace
namespace qpixl {
  
  /// Return whether x is zero or a power of 2
  static inline bool ispow2( const size_t x ) {
      return  !(x & (x-1));
  }

  /// Returns the exponent of a power of 2.
  static inline size_t ilog2( const size_t x ) {
    assert( std::popcount( x ) == 1 ) ;
    return std::countr_zero( x ) ;
  }

  /// Return the next power of two
  static inline size_t nextpow2( size_t x ) {
    x-- ;
    x |= x >>  1 ;
    x |= x >>  2 ;
    x |= x >>  4 ;
    x |= x >>  8 ;
    x |= x >> 16 ;
    x |= x >> 32 ;
    x++ ;
    return x ;
  }

}

#endif
