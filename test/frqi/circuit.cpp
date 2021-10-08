#include <gtest/gtest.h>
#include "qclab/QCircuit.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "qpixl/frqi/util.hpp"
#include "qpixl/frqi/circuit.hpp"

template < typename MAT1, typename MAT2, typename R >
void test_qpixl_frqi_check_matrix( const MAT1& m1, const MAT2& m2, const R eps ) {
  assert( m1.size() == m2.size() ) ;
  for ( size_t i = 0; i < m1.size(); i++ ) {
    for ( size_t j = 0; j < m1.size(); j++ ) {
      EXPECT_NEAR( m1(i,j) , m2(i,j)  , eps );
    }
  }
}

template <typename R>
void test_qpixl_frqi_compressedFRQICircuit() {
  const R eps = std::numeric_limits< R >::epsilon() ;

  std::vector< R > a = { 0.5 , 1 } ;

  auto circA = qpixl::frqi::compressedFRQICircuit( a ) ;

  // matrix
  qclab::dense::SquareMatrix< R > circA_check(    
  0.620544580563746, -0.339005049421045,  0.620544580563746, -0.339005049421045,
  0.339005049421045,  0.620544580563746,  0.339005049421045,  0.620544580563746,
  0.382051424370090, -0.595009839529386, -0.382051424370090,  0.595009839529386,
  0.595009839529386,  0.382051424370090, -0.595009839529386, -0.382051424370090
  ) ;

  test_qpixl_frqi_check_matrix( circA.matrix() , circA_check  , 10*eps );
  
  
  std::vector< R > b = { 0.5 , 1 , 1.5 , 2 };

  auto circB = qpixl::frqi::compressedFRQICircuit( b ) ;

  // matrix
  qclab::dense::SquareMatrix< R > circB_check( 8 ) ;

  circB_check( 0, 0 ) =  0.438791280945186 ;
  circB_check( 1, 0 ) =  0.239712769302101 ;
  circB_check( 2, 0 ) =  0.270151152934070 ;
  circB_check( 3, 0 ) =  0.420735492403948 ;
  circB_check( 4, 0 ) =  0.035368600833851 ;
  circB_check( 5, 0 ) =  0.498747493302027 ;
  circB_check( 6, 0 ) = -0.208073418273571 ;
  circB_check( 7, 0 ) =  0.454648713412841 ;

  circB_check( 0, 1 ) = -0.239712769302101 ;
  circB_check( 1, 1 ) =  0.438791280945186 ;
  circB_check( 2, 1 ) = -0.420735492403948 ;
  circB_check( 3, 1 ) =  0.270151152934070 ;
  circB_check( 4, 1 ) = -0.498747493302027 ;
  circB_check( 5, 1 ) =  0.035368600833851 ;
  circB_check( 6, 1 ) = -0.454648713412841 ;
  circB_check( 7, 1 ) = -0.208073418273571 ;

  circB_check( 0, 2 ) =  0.438791280945186 ;
  circB_check( 1, 2 ) =  0.239712769302101 ;
  circB_check( 2, 2 ) = -0.270151152934070 ;
  circB_check( 3, 2 ) = -0.420735492403948 ;
  circB_check( 4, 2 ) =  0.035368600833851 ;
  circB_check( 5, 2 ) =  0.498747493302027 ;
  circB_check( 6, 2 ) =  0.208073418273571 ;
  circB_check( 7, 2 ) = -0.454648713412841 ;

  circB_check( 0, 3 ) = -0.239712769302101 ;
  circB_check( 1, 3 ) =  0.438791280945186 ;
  circB_check( 2, 3 ) =  0.420735492403948 ;
  circB_check( 3, 3 ) = -0.270151152934070 ;
  circB_check( 4, 3 ) = -0.498747493302027 ;
  circB_check( 5, 3 ) =  0.035368600833851 ;
  circB_check( 6, 3 ) =  0.454648713412841 ;
  circB_check( 7, 3 ) =  0.208073418273571 ;

  circB_check( 0, 4 ) =  0.438791280945186 ;
  circB_check( 1, 4 ) =  0.239712769302101 ;
  circB_check( 2, 4 ) =  0.270151152934070 ;
  circB_check( 3, 4 ) =  0.420735492403948 ;
  circB_check( 4, 4 ) = -0.035368600833851 ;
  circB_check( 5, 4 ) = -0.498747493302027 ;
  circB_check( 6, 4 ) =  0.208073418273571 ;
  circB_check( 7, 4 ) = -0.454648713412841 ;

  circB_check( 0, 5 ) = -0.239712769302101 ;
  circB_check( 1, 5 ) =  0.438791280945186 ;
  circB_check( 2, 5 ) = -0.420735492403948 ;
  circB_check( 3, 5 ) =  0.270151152934070 ;
  circB_check( 4, 5 ) =  0.498747493302027 ;
  circB_check( 5, 5 ) = -0.035368600833851 ;
  circB_check( 6, 5 ) =  0.454648713412841 ;
  circB_check( 7, 5 ) =  0.208073418273571 ;

  circB_check( 0, 6 ) =  0.438791280945186 ;
  circB_check( 1, 6 ) =  0.239712769302101 ;
  circB_check( 2, 6 ) = -0.270151152934070 ;
  circB_check( 3, 6 ) = -0.420735492403948 ;
  circB_check( 4, 6 ) = -0.035368600833851 ;
  circB_check( 5, 6 ) = -0.498747493302027 ;
  circB_check( 6, 6 ) = -0.208073418273571 ;
  circB_check( 7, 6 ) =  0.454648713412841 ;

  circB_check( 0, 7 ) = -0.239712769302101 ;
  circB_check( 1, 7 ) =  0.438791280945186 ;
  circB_check( 2, 7 ) =  0.420735492403948 ;
  circB_check( 3, 7 ) = -0.270151152934070 ;
  circB_check( 4, 7 ) =  0.498747493302027 ;
  circB_check( 5, 7 ) = -0.035368600833851 ;
  circB_check( 6, 7 ) = -0.454648713412841 ;
  circB_check( 7, 7 ) = -0.208073418273571 ;


  test_qpixl_frqi_check_matrix( circB.matrix() , circB_check  , 10*eps );



  std::vector< R > c = { 0.5 , 1 , 1.5 , 2 , 2.5 , 3 , 3.5 , 4 };

  auto circC = qpixl::frqi::compressedFRQICircuit( c ) ;  

  // matrix
  qclab::dense::SquareMatrix< R > circC_check( 16 ) ;

  circC_check(  0, 0 ) =    0.310272290281873 ;
  circC_check(  1, 0 ) =    0.169502524710522 ;
  circC_check(  2, 0 ) =    0.191025712185045 ;
  circC_check(  3, 0 ) =    0.297504919764693 ;
  circC_check(  4, 0 ) =    0.025009377490697 ;
  circC_check(  5, 0 ) =    0.352667734613656 ;
  circC_check(  6, 0 ) =   -0.147130125045907 ;
  circC_check(  7, 0 ) =    0.321485188311959 ;
  circC_check(  8, 0 ) =   -0.283247041628773 ;
  circC_check(  9, 0 ) =    0.211591855723580 ;
  circC_check( 10, 0 ) =   -0.350015203834987 ;
  circC_check( 11, 0 ) =    0.049893457330116 ;
  circC_check( 12, 0 ) =   -0.331087436935406 ;
  circC_check( 13, 0 ) =   -0.124020599512917 ;
  circC_check( 14, 0 ) =   -0.231097918395994 ;
  circC_check( 15, 0 ) =   -0.267570088225568 ;

  circC_check(  0, 1 ) =   -0.169502524710522 ;
  circC_check(  1, 1 ) =    0.310272290281873 ;
  circC_check(  2, 1 ) =   -0.297504919764693 ;
  circC_check(  3, 1 ) =    0.191025712185045 ;
  circC_check(  4, 1 ) =   -0.352667734613656 ;
  circC_check(  5, 1 ) =    0.025009377490697 ;
  circC_check(  6, 1 ) =   -0.321485188311959 ;
  circC_check(  7, 1 ) =   -0.147130125045907 ;
  circC_check(  8, 1 ) =   -0.211591855723580 ;
  circC_check(  9, 1 ) =   -0.283247041628773 ;
  circC_check( 10, 1 ) =   -0.049893457330116 ;
  circC_check( 11, 1 ) =   -0.350015203834987 ;
  circC_check( 12, 1 ) =    0.124020599512917 ;
  circC_check( 13, 1 ) =   -0.331087436935406 ;
  circC_check( 14, 1 ) =    0.267570088225568 ;
  circC_check( 15, 1 ) =   -0.231097918395994 ;

  circC_check(  0, 2 ) =    0.310272290281873 ;
  circC_check(  1, 2 ) =    0.169502524710522 ;
  circC_check(  2, 2 ) =   -0.191025712185045 ;
  circC_check(  3, 2 ) =   -0.297504919764693 ;
  circC_check(  4, 2 ) =    0.025009377490697 ;
  circC_check(  5, 2 ) =    0.352667734613656 ;
  circC_check(  6, 2 ) =    0.147130125045907 ;
  circC_check(  7, 2 ) =   -0.321485188311959 ;
  circC_check(  8, 2 ) =   -0.283247041628773 ;
  circC_check(  9, 2 ) =    0.211591855723580 ;
  circC_check( 10, 2 ) =    0.350015203834987 ;
  circC_check( 11, 2 ) =   -0.049893457330116 ;
  circC_check( 12, 2 ) =   -0.331087436935406 ;
  circC_check( 13, 2 ) =   -0.124020599512917 ;
  circC_check( 14, 2 ) =    0.231097918395994 ;
  circC_check( 15, 2 ) =    0.267570088225568 ;

  circC_check(  0, 3 ) =   -0.169502524710522 ;
  circC_check(  1, 3 ) =    0.310272290281873 ;
  circC_check(  2, 3 ) =    0.297504919764693 ;
  circC_check(  3, 3 ) =   -0.191025712185045 ;
  circC_check(  4, 3 ) =   -0.352667734613656 ;
  circC_check(  5, 3 ) =    0.025009377490697 ;
  circC_check(  6, 3 ) =    0.321485188311959 ;
  circC_check(  7, 3 ) =    0.147130125045907 ;
  circC_check(  8, 3 ) =   -0.211591855723580 ;
  circC_check(  9, 3 ) =   -0.283247041628773 ;
  circC_check( 10, 3 ) =    0.049893457330116 ;
  circC_check( 11, 3 ) =    0.350015203834987 ;
  circC_check( 12, 3 ) =    0.124020599512917 ;
  circC_check( 13, 3 ) =   -0.331087436935406 ;
  circC_check( 14, 3 ) =   -0.267570088225568 ;
  circC_check( 15, 3 ) =    0.231097918395994 ;

  circC_check(  0, 4 ) =    0.310272290281873 ;
  circC_check(  1, 4 ) =    0.169502524710522 ;
  circC_check(  2, 4 ) =    0.191025712185045 ;
  circC_check(  3, 4 ) =    0.297504919764693 ;
  circC_check(  4, 4 ) =   -0.025009377490697 ;
  circC_check(  5, 4 ) =   -0.352667734613656 ;
  circC_check(  6, 4 ) =    0.147130125045907 ;
  circC_check(  7, 4 ) =   -0.321485188311959 ;
  circC_check(  8, 4 ) =   -0.283247041628773 ;
  circC_check(  9, 4 ) =    0.211591855723580 ;
  circC_check( 10, 4 ) =   -0.350015203834987 ;
  circC_check( 11, 4 ) =    0.049893457330116 ;
  circC_check( 12, 4 ) =    0.331087436935406 ;
  circC_check( 13, 4 ) =    0.124020599512917 ;
  circC_check( 14, 4 ) =    0.231097918395994 ;
  circC_check( 15, 4 ) =    0.267570088225568 ;

  circC_check(  0, 5 ) =   -0.169502524710522 ;
  circC_check(  1, 5 ) =    0.310272290281873 ;
  circC_check(  2, 5 ) =   -0.297504919764693 ;
  circC_check(  3, 5 ) =    0.191025712185045 ;
  circC_check(  4, 5 ) =    0.352667734613656 ;
  circC_check(  5, 5 ) =   -0.025009377490697 ;
  circC_check(  6, 5 ) =    0.321485188311959 ;
  circC_check(  7, 5 ) =    0.147130125045907 ;
  circC_check(  8, 5 ) =   -0.211591855723580 ;
  circC_check(  9, 5 ) =   -0.283247041628773 ;
  circC_check( 10, 5 ) =   -0.049893457330116 ;
  circC_check( 11, 5 ) =   -0.350015203834987 ;
  circC_check( 12, 5 ) =   -0.124020599512917 ;
  circC_check( 13, 5 ) =    0.331087436935406 ;
  circC_check( 14, 5 ) =   -0.267570088225568 ;
  circC_check( 15, 5 ) =    0.231097918395994 ;

  circC_check(  0, 6 ) =    0.310272290281873 ;
  circC_check(  1, 6 ) =    0.169502524710522 ;
  circC_check(  2, 6 ) =   -0.191025712185045 ;
  circC_check(  3, 6 ) =   -0.297504919764693 ;
  circC_check(  4, 6 ) =   -0.025009377490697 ;
  circC_check(  5, 6 ) =   -0.352667734613656 ;
  circC_check(  6, 6 ) =   -0.147130125045907 ;
  circC_check(  7, 6 ) =    0.321485188311959 ;
  circC_check(  8, 6 ) =   -0.283247041628773 ;
  circC_check(  9, 6 ) =    0.211591855723580 ;
  circC_check( 10, 6 ) =    0.350015203834987 ;
  circC_check( 11, 6 ) =   -0.049893457330116 ;
  circC_check( 12, 6 ) =    0.331087436935406 ;
  circC_check( 13, 6 ) =    0.124020599512917 ;
  circC_check( 14, 6 ) =   -0.231097918395994 ;
  circC_check( 15, 6 ) =   -0.267570088225568 ;

  circC_check(  0, 7 ) =   -0.169502524710522 ;
  circC_check(  1, 7 ) =    0.310272290281873 ;
  circC_check(  2, 7 ) =    0.297504919764693 ;
  circC_check(  3, 7 ) =   -0.191025712185045 ;
  circC_check(  4, 7 ) =    0.352667734613656 ;
  circC_check(  5, 7 ) =   -0.025009377490697 ;
  circC_check(  6, 7 ) =   -0.321485188311959 ;
  circC_check(  7, 7 ) =   -0.147130125045907 ;
  circC_check(  8, 7 ) =   -0.211591855723580 ;
  circC_check(  9, 7 ) =   -0.283247041628773 ;
  circC_check( 10, 7 ) =    0.049893457330116 ;
  circC_check( 11, 7 ) =    0.350015203834987 ;
  circC_check( 12, 7 ) =   -0.124020599512917 ;
  circC_check( 13, 7 ) =    0.331087436935406 ;
  circC_check( 14, 7 ) =    0.267570088225568 ;
  circC_check( 15, 7 ) =   -0.231097918395994 ;

  circC_check(  0, 8 ) =    0.310272290281873 ;
  circC_check(  1, 8 ) =    0.169502524710522 ;
  circC_check(  2, 8 ) =    0.191025712185045 ;
  circC_check(  3, 8 ) =    0.297504919764693 ;
  circC_check(  4, 8 ) =    0.025009377490697 ;
  circC_check(  5, 8 ) =    0.352667734613656 ;
  circC_check(  6, 8 ) =   -0.147130125045907 ;
  circC_check(  7, 8 ) =    0.321485188311959 ;
  circC_check(  8, 8 ) =    0.283247041628773 ;
  circC_check(  9, 8 ) =   -0.211591855723580 ;
  circC_check( 10, 8 ) =    0.350015203834987 ;
  circC_check( 11, 8 ) =   -0.049893457330116 ;
  circC_check( 12, 8 ) =    0.331087436935406 ;
  circC_check( 13, 8 ) =    0.124020599512917 ;
  circC_check( 14, 8 ) =    0.231097918395994 ;
  circC_check( 15, 8 ) =    0.267570088225568 ;

  circC_check(  0, 9 ) =   -0.169502524710522 ;
  circC_check(  1, 9 ) =    0.310272290281873 ;
  circC_check(  2, 9 ) =   -0.297504919764693 ;
  circC_check(  3, 9 ) =    0.191025712185045 ;
  circC_check(  4, 9 ) =   -0.352667734613656 ;
  circC_check(  5, 9 ) =    0.025009377490697 ;
  circC_check(  6, 9 ) =   -0.321485188311959 ;
  circC_check(  7, 9 ) =   -0.147130125045907 ;
  circC_check(  8, 9 ) =    0.211591855723580 ;
  circC_check(  9, 9 ) =    0.283247041628773 ;
  circC_check( 10, 9 ) =    0.049893457330116 ;
  circC_check( 11, 9 ) =    0.350015203834987 ;
  circC_check( 12, 9 ) =   -0.124020599512917 ;
  circC_check( 13, 9 ) =    0.331087436935406 ;
  circC_check( 14, 9 ) =   -0.267570088225568 ;
  circC_check( 15, 9 ) =    0.231097918395994 ;

  circC_check(  0, 10 ) =    0.310272290281873 ;
  circC_check(  1, 10 ) =    0.169502524710522 ;
  circC_check(  2, 10 ) =   -0.191025712185045 ;
  circC_check(  3, 10 ) =   -0.297504919764693 ;
  circC_check(  4, 10 ) =    0.025009377490697 ;
  circC_check(  5, 10 ) =    0.352667734613656 ;
  circC_check(  6, 10 ) =    0.147130125045907 ;
  circC_check(  7, 10 ) =   -0.321485188311959 ;
  circC_check(  8, 10 ) =    0.283247041628773 ;
  circC_check(  9, 10 ) =   -0.211591855723580 ;
  circC_check( 10, 10 ) =   -0.350015203834987 ;
  circC_check( 11, 10 ) =    0.049893457330116 ;
  circC_check( 12, 10 ) =    0.331087436935406 ;
  circC_check( 13, 10 ) =    0.124020599512917 ;
  circC_check( 14, 10 ) =   -0.231097918395994 ;
  circC_check( 15, 10 ) =   -0.267570088225568 ;

  circC_check(  0, 11 ) =   -0.169502524710522 ;
  circC_check(  1, 11 ) =    0.310272290281873 ;
  circC_check(  2, 11 ) =    0.297504919764693 ;
  circC_check(  3, 11 ) =   -0.191025712185045 ;
  circC_check(  4, 11 ) =   -0.352667734613656 ;
  circC_check(  5, 11 ) =    0.025009377490697 ;
  circC_check(  6, 11 ) =    0.321485188311959 ;
  circC_check(  7, 11 ) =    0.147130125045907 ;
  circC_check(  8, 11 ) =    0.211591855723580 ;
  circC_check(  9, 11 ) =    0.283247041628773 ;
  circC_check( 10, 11 ) =   -0.049893457330116 ;
  circC_check( 11, 11 ) =   -0.350015203834987 ;
  circC_check( 12, 11 ) =   -0.124020599512917 ;
  circC_check( 13, 11 ) =    0.331087436935406 ;
  circC_check( 14, 11 ) =    0.267570088225568 ;
  circC_check( 15, 11 ) =   -0.231097918395994 ;

  circC_check(  0, 12 ) =    0.310272290281873 ;
  circC_check(  1, 12 ) =    0.169502524710522 ;
  circC_check(  2, 12 ) =    0.191025712185045 ;
  circC_check(  3, 12 ) =    0.297504919764693 ;
  circC_check(  4, 12 ) =   -0.025009377490697 ;
  circC_check(  5, 12 ) =   -0.352667734613656 ;
  circC_check(  6, 12 ) =    0.147130125045907 ;
  circC_check(  7, 12 ) =   -0.321485188311959 ;
  circC_check(  8, 12 ) =    0.283247041628773 ;
  circC_check(  9, 12 ) =   -0.211591855723580 ;
  circC_check( 10, 12 ) =    0.350015203834987 ;
  circC_check( 11, 12 ) =   -0.049893457330116 ;
  circC_check( 12, 12 ) =   -0.331087436935406 ;
  circC_check( 13, 12 ) =   -0.124020599512917 ;
  circC_check( 14, 12 ) =   -0.231097918395994 ;
  circC_check( 15, 12 ) =   -0.267570088225568 ;

  circC_check(  0, 13 ) =   -0.169502524710522 ;
  circC_check(  1, 13 ) =    0.310272290281873 ;
  circC_check(  2, 13 ) =   -0.297504919764693 ;
  circC_check(  3, 13 ) =    0.191025712185045 ;
  circC_check(  4, 13 ) =    0.352667734613656 ;
  circC_check(  5, 13 ) =   -0.025009377490697 ;
  circC_check(  6, 13 ) =    0.321485188311959 ;
  circC_check(  7, 13 ) =    0.147130125045907 ;
  circC_check(  8, 13 ) =    0.211591855723580 ;
  circC_check(  9, 13 ) =    0.283247041628773 ;
  circC_check( 10, 13 ) =    0.049893457330116 ;
  circC_check( 11, 13 ) =    0.350015203834987 ;
  circC_check( 12, 13 ) =    0.124020599512917 ;
  circC_check( 13, 13 ) =   -0.331087436935406 ;
  circC_check( 14, 13 ) =    0.267570088225568 ;
  circC_check( 15, 13 ) =   -0.231097918395994 ;

  circC_check(  0, 14 ) =    0.310272290281873 ;
  circC_check(  1, 14 ) =    0.169502524710522 ;
  circC_check(  2, 14 ) =   -0.191025712185045 ;
  circC_check(  3, 14 ) =   -0.297504919764693 ;
  circC_check(  4, 14 ) =   -0.025009377490697 ;
  circC_check(  5, 14 ) =   -0.352667734613656 ;
  circC_check(  6, 14 ) =   -0.147130125045907 ;
  circC_check(  7, 14 ) =    0.321485188311959 ;
  circC_check(  8, 14 ) =    0.283247041628773 ;
  circC_check(  9, 14 ) =   -0.211591855723580 ;
  circC_check( 10, 14 ) =   -0.350015203834987 ;
  circC_check( 11, 14 ) =    0.049893457330116 ;
  circC_check( 12, 14 ) =   -0.331087436935406 ;
  circC_check( 13, 14 ) =   -0.124020599512917 ;
  circC_check( 14, 14 ) =    0.231097918395994 ;
  circC_check( 15, 14 ) =    0.267570088225568 ;

  circC_check(  0, 15 ) =   -0.169502524710522 ;
  circC_check(  1, 15 ) =    0.310272290281873 ;
  circC_check(  2, 15 ) =    0.297504919764693 ;
  circC_check(  3, 15 ) =   -0.191025712185045 ;
  circC_check(  4, 15 ) =    0.352667734613656 ;
  circC_check(  5, 15 ) =   -0.025009377490697 ;
  circC_check(  6, 15 ) =   -0.321485188311959 ;
  circC_check(  7, 15 ) =   -0.147130125045907 ;
  circC_check(  8, 15 ) =    0.211591855723580 ;
  circC_check(  9, 15 ) =    0.283247041628773 ;
  circC_check( 10, 15 ) =   -0.049893457330116 ;
  circC_check( 11, 15 ) =   -0.350015203834987 ;
  circC_check( 12, 15 ) =    0.124020599512917 ;
  circC_check( 13, 15 ) =   -0.331087436935406 ;
  circC_check( 14, 15 ) =   -0.267570088225568 ;
  circC_check( 15, 15 ) =    0.231097918395994 ;

  test_qpixl_frqi_check_matrix( circC.matrix() , circC_check  , 10*eps );
}

template <typename R>
void test_qpixl_frqi_circuit() {
  test_qpixl_frqi_compressedFRQICircuit< R >() ;
}

/*
 * float
 */

TEST( qpixl_frqi_circuit , float ) {
  test_qpixl_frqi_circuit< float >() ;
}


/*
 * double
 */
TEST( qpixl_frqi_circuit , double ) {
  test_qpixl_frqi_circuit< double >() ;
}