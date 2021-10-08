// (C) Copyright Daan Camps, Mercy Amankwah, E. Wes Bethel, Talita Perciano and
//               Roel Van Beeumen

#ifndef qpixl_frqi_circuit_hpp
#define qpixl_frqi_circuit_hpp

#include "qclab/QCircuit.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "qpixl/util.hpp"
#include "qpixl/frqi/util.hpp"
#include <cmath>

/// QPIXL::FRQI namespace
namespace qpixl::frqi {

  /**
   * \brief Generates a compressed FRQI circuit for the input vector `a`
   * 
   *   - `a` contains the angle representation of the image, its size should
   *     be a power of 2.
   * 
   *   - `compression` contains the compression value between 0 and 100 with 0
   *     meaning no compression at all ( all image coefficients are encoded ) 
   *     and 100 meaning full compression ( no image coefficients are encoded )
   *
   */
  template <typename V>
  auto compressedFRQICircuit( V& a , const double compression = 0.0 ) {
    using T = typename V::value_type ;
    using H  = qclab::qgates::Hadamard< T > ;
    using CNOT = qclab::qgates::CNOT< T > ;
    using RY = qclab::qgates::RotationY< T > ;

    size_t n = a.size() ;
    assert( qpixl::ispow2( n ) ) ;
    size_t k = qpixl::ilog2( n ) ;

    // Multiply angles by two
    std::for_each( a.begin(), a.end(), []( T& n ){ n *= 2; } ) ;

    // Convert angles through a scaled permuted fast Walsh-Hadamard transform
    qpixl::frqi::sfwht( a ) ;
    qpixl::frqi::grayPermutation( a ) ;
    
    // idx = [0,1,2,...,n-1]                                             
    std::vector< size_t > index( n ) ;                        
    std::iota( index.begin() , index.end() , 0 ) ;                              
    
    // sort vector a by absolute values and get ordering in index array
    std::sort( index.begin() , index.end() ,                                    
               [&a]( size_t i1 , size_t i2 ) 
               { return std::abs(a[i1]) < std::abs(a[i2]); } ) ;

    // set smallest absolute values of a to zero according to compression param
    size_t cutoff = size_t( (compression/100.0) * n );
    for ( std::vector<size_t>::const_iterator it = index.begin(); 
          it < index.begin() + cutoff; ++it ) {
      a[*it] = 0.0 ;
    }    

    // Construct FRQI circuit
    qclab::QCircuit< T > circuit( k + 1 ) ;
    
    // Hadamard register
    for ( size_t i = 0; i < k; i++ ) {
      circuit.push_back( std::make_unique< H >( i ) ) ;
    }

    // Compressed uniformly controlled rotation register
    size_t ctrl, pc, i = 0 ;
    while ( i < (size_t(1)<<k) ) {
      
      // Reset parity check
      pc = 0 ;

      // Add RY gate
      if ( a[i] != 0 ) {
        circuit.push_back( std::make_unique< RY >( int(k) , a[i] ) ) ;
      }

      // Loop over sequence of consecutive zero angles
      do {

        // Compute control qubit
        if ( i ==  (size_t(1)<<k) - 1 ) {
          ctrl = 0 ;
        } else {
          ctrl = qpixl::frqi::grayCode( i ) ^ qpixl::frqi::grayCode( i + 1 ) ;
          assert( std::popcount( ctrl ) == 1 ) ;
          ctrl = k - std::countr_zero( ctrl ) - 1 ;  
        }

        // Update parity check
        pc ^= 1UL << ctrl ;

        i++ ;
      } while ( i < (size_t(1)<<k) && a[i] == 0 ) ;

      // Add CNOTs determined by parity check
      for ( int j = 0; j < k; j++ ) {
        if ( ( pc >> j ) & 1U ) {          
          circuit.push_back( std::make_unique< CNOT >( j , k ) ) ;
        }
      }      
    }
        
    // return output
    return circuit ;
  }

} // namespace qpixl::frqi

#endif