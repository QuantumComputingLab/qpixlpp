// (C) Copyright Daan Camps, Mercy Amankwah, E. Wes Bethel, Talita Perciano and
//               Roel Van Beeumen

#include "qpixl/frqi/util.hpp"
#include <chrono>
#include <random>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using TP = std::chrono::time_point< std::chrono::high_resolution_clock > ;

void tic( TP& t_bgn ) {
  t_bgn = std::chrono::high_resolution_clock::now() ;
} // tic

void tic( const std::string name , TP& t_bgn ) {
  std::cout << name << std::endl ;
  t_bgn = std::chrono::high_resolution_clock::now() ;
} // tic

double toc( const TP& t_bgn ) {
  TP t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  return time ;
} // toc

double toc( const TP& t_bgn , TP& t_end ) {
  t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  std::cout << "                                   " << time << "s\n" ;
  return time ;
} // toc

template <typename T>
auto pow2Vector( const size_t k ) {

  // random entries
  std::random_device                rd ;
  std::mt19937_64                   gen( rd() ) ;
  std::uniform_real_distribution<>  dis( 0.0, 1.0 ) ;

  // generate vector
  std::vector< T > v( size_t(1)<<k ) ;
  for ( size_t i = 0; i < (size_t(1)<<k); i++ ) {
    v[i] = dis(gen) ;
  }

  return v ;
}

template <typename T>
int timings( const std::vector< size_t>& N, const int outer , const int inner ){
  std::cout << "outer = " << outer << ", inner = " << inner << std::endl ;
  //problem dimension
  std::cout << "qubits = " << N.front() << ":" << N.back() ;
  #ifdef _OPENMP
  std::cout << ", omp_get_max_threads() = " << omp_get_max_threads() ;
  #endif
  std::cout << std::endl << std::endl ;

  // timing variables
  TP  t_bgn ;
  TP  t_end ;

  std::vector< int >     qubits ;
  std::vector< double >  time_sfwht, time_gray ;
  for ( size_t n : N ) {

    std::cout << "n = " << n << ":" << std::endl ;
    qubits.push_back( n ) ;
    
    // outer loop
    for ( int o = 0; o < outer; o++ ) {
      
      // generate random vectors
      tic( "  * Constructing random vectors..." , t_bgn ) ;
      std::vector< std::vector< T > > vectors ;
      for ( int i = 0; i < inner; i++ ) {
        vectors.push_back( pow2Vector< T >( n ) );
      }
      toc( t_bgn , t_end ) ;

      // scaled fast Walsh-Hadamard transform
      tic( "  * fast Walsh-Hadamard transform..." , t_bgn ) ;
      for ( int i = 0; i < inner; i++ ) {
        qpixl::frqi::sfwht( vectors[i] );
      }
      const double ttot_sfwht = toc( t_bgn , t_end ) / inner ;
      if ( o == 0 ) time_sfwht.push_back( ttot_sfwht ) ;
      else if ( ttot_sfwht < time_sfwht.back() ) time_sfwht.back() = 
        ttot_sfwht ;

      // in-place Gray permutation
      tic( "  * Gray permutation..." , t_bgn ) ;
      for ( int i = 0; i < inner; i++ ) {
        qpixl::frqi::grayPermutation( vectors[i] );
      }
      const double ttot_gray = toc( t_bgn , t_end ) / inner ;
      if ( o == 0 ) time_gray.push_back( ttot_gray ) ;
      else if ( ttot_gray < time_gray.back() ) time_gray.back() = ttot_gray ;
    }
  }

  // output
  std::cout << "results = [" << std::endl ;
  for ( size_t i = 0; i < time_sfwht.size(); i++ ) {
    std::printf( "%6i, %10.4e, %10.4e" , qubits[i] , time_sfwht[i], 
      time_gray[i] ) ;
    if ( i == time_sfwht.size() - 1 ) {
      std::printf( "];\n" ) ;
    } else {
      std::printf( " ;\n" ) ;
    }
  }
  std::cout << std::endl ;

  // successful
  return 0 ;
}

int main( int argc , char *argv[] ) {

  // parse user input
  char type  = 'd' ; if ( argc > 1 ) type  = argv[1][0] ;
  int  nmin  = 3 ;   if ( argc > 2 ) nmin  = std::stoi( argv[2] ) ;
  int  nmax  = 20 ;  if ( argc > 3 ) nmax  = std::stoi( argv[3] ) ;
  int  step  = 1 ;   if ( argc > 4 ) step  = std::stoi( argv[4] ) ;
  int  outer = 1 ;   if ( argc > 5 ) outer = std::stoi( argv[5] ) ;
  int  inner = 1 ;   if ( argc > 6 ) inner = std::stoi( argv[6] ) ;

  std::vector< size_t > N ;
  // linear scale
  for ( int i = nmin; i <= nmax; i += step ) {
    N.push_back( i ) ;
  }
  if ( N.back() != nmax ) { N.push_back( nmax ) ; }

  int r = 0 ;
  if ( type == 's' ) {
    // float
    std::cout << "\n+++ SFWHT and Gray permutation (float) +++\n\n" ;
    r = timings< float >( N , outer , inner ) ;
    if ( r != 0 ) return r ;
  } else if ( type == 'd' ) {
    // double
    std::cout << "\n+++ SFWHT and Gray permutation (double) +++\n\n" ;
    r = timings< double >( N , outer , inner ) ;
    if ( r != 0 ) return r + 1000 ;
  } else {
    return -2 ;
  }

  // successful
  std::cout << ">> end SFWHT and Gray permutation <<" << std::endl ;

  return 0 ;

}