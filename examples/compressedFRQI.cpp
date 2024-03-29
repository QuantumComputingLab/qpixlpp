// (C) Copyright Daan Camps, Mercy Amankwah, E. Wes Bethel, Talita Perciano and
//               Roel Van Beeumen

#include "qpixl/frqi/util.hpp"
#include "qpixl/util.hpp"
#include "qpixl/pgm.hpp"
#include "qpixl/frqi/circuit.hpp"
#include "qclab/QCircuit.hpp"
#include <cmath>
#include <iomanip>

int main( int argc , char *argv[] ) {

  // parse input arguments
  if ( argc == 2 ) {
    if ( std::string(argv[1]) == "--help" ) {
      std::cout << std::endl
                << "QPIXL++: Compile quantum circuits for image representations"
                << std::endl
                << "***********************************************************"
                << std::endl << std::endl
                << "compressedFRQI requires two mandatory input arguments, "
                << std::endl << "and two optional input arguments."
                << std::endl << std::endl
                << "1. path to input PGM image file. Image data is expected"
                << std::endl
                << "   in ASCII portable gray map format without comment."
                << std::endl << std::endl
                << "2. path and name structure for all output data."
                << std::endl
                << "   For example : " << std::endl
                << "       /mypath/output " << std::endl
                << "       This will generate output.qasm, output.pgm, and "
                << std::endl
                << "       output_sim.pgm (depending on optional parameters)"
                << std::endl
                << "       in /mypath/"
                << std::endl << std::endl
                << "3. (optional) compression level between 0 and 100 percent"
                << std::endl
                << "   default compression level is 0."
                << std::endl << std::endl
                << "4. (optional) simulate flag: 0 or 1 to determine if circuit"
                << std::endl
                << "   is simulated with QCLAB++, default is 0."
                << std::endl << std::endl ;
      return 0 ;
    }
  }
  if ( argc < 3 ) {
    std::cerr << "ERROR -- input and output filenames are required!" 
              << std::endl 
              << "Use --help for more more ."
              << std::endl ;
    return -1 ;
  }
  std::string fileNameIn(  argv[1] ) ;
  std::string fileNameOut( argv[2] ) ;
  double  compr  = 0 ; if ( argc > 3 ) compr  = std::stod( argv[3] ) ;
  int sim = 0 ; if ( argc > 4 ) sim = std::stoi( argv[4] ) ;

  // Read grayscale image data
  std::vector< double > data ;
  size_t nrows, ncols, maxval ;
  int s, zeroPad = 2 ;
  std::cout << " * Reading image data from file.. " << std::endl ;
  s = qpixl::readPGMA(fileNameIn , data , nrows , ncols , maxval , zeroPad) ;

  if ( s != 0 ) {
    std::cerr << "ERROR -- reading input image file" << std::endl;
    return -2 ;
  }

  // Image statistics
  std::cout << std::endl
              << "   Image statistics : " << std::endl
              << "     * number of rows             : " << nrows << std::endl
              << "     * number of columns          : " << ncols << std::endl
              << "     * maximum pixel value        : " << maxval << std::endl
              << "     * size of zero padded vector : " << data.size() 
              << std::endl << std::endl ;
    
  // Convert data to angles for FRQI
  std::cout << " * Converting grayscale to angles.. " << std::endl ;
  qpixl::frqi::convertToAngles( data , maxval ) ;

  // Create compressed FRQI circuit
  std::cout << " * Compressing FRQI circuit.. " << std::endl ;
  auto circuit = qpixl::frqi::compressedFRQICircuit( data , compr ) ;

  // Count CNOTs and single qubit gates in FRQI circuit
  int nCNOT = 0, nRY = 0 ;
  const int nH = circuit.nbQubits() - 1 ;
  using citer  = typename qclab::QCircuit< double >::const_iterator ;
  for ( citer it = circuit.begin(); it !=circuit.end(); it++ ) {
    if ( (*it)->nbQubits() == 1 ) nRY += 1 ;
    else nCNOT += 1 ;
  }
  nRY -= nH ;
  
  const double maxgates = int(1)<<(circuit.nbQubits()-1) ;

  // Write QASM to file
  std::cout << " * Writing compressed FRQI circuit to QASM.. " << std::endl ;
  std::stringstream qasm ;
  circuit.toQASM( qasm ) ;
  std::string fileNameQASM ;
  fileNameQASM.append( fileNameOut ) ;
  fileNameQASM.append( ".qasm" ) ;
  std::ofstream ostrm( fileNameQASM , std::ios::out ) ;
  const int N = circuit.nbQubits() ;
  ostrm << "// Generated by QPIXL++" << std::endl
        << "// https://github.com/QuantumComputingLab/qpixlpp" << std::endl
        << std::endl << std::endl
        << "//   Circuit statistics : " << std::endl
        << "//     * number of qubits               : " 
        << circuit.nbQubits() << std::endl
        << "//     * total number of gates          : " << circuit.nbGates()
        << std::endl
        << "//       - number of CNOTs              : " << nCNOT 
        << std::endl
        << "//       - number of Ry gates           : " << nRY 
        << std::endl
        << "//       - number of Hadamard gates     : " << nH 
        << std::endl
        << "//     * Compression setting            : " << compr
        << std::endl
        << "//       - CNOT compression ratio [%]   : " 
        << (1.0 - (double(nCNOT)/maxgates))*100
        << std::endl
        << "//       - Ry compression ratio [%]     : " 
        << (1.0 - (double(nRY)/maxgates))*100
        << std::endl << std::endl
        << "OPENQASM 2.0;" << std::endl
        << "include \"qelib1.inc\";" << std::endl << std::endl
        << "qreg q[" << N << "];" << std::endl
        << qasm.str() ;
  ostrm.close() ;

  // Circuit statistics
  std::cout << std::fixed ;
  std::cout << std::setprecision(1) ;
  std::cout << std::endl
            << "   Circuit statistics : " << std::endl
            << "     * number of qubits               : " 
            << circuit.nbQubits() << std::endl
            << "     * total number of gates          : " << circuit.nbGates()
            << std::endl
            << "       - number of CNOTs              : " << nCNOT 
            << std::endl
            << "       - number of Ry gates           : " << nRY
            << std::endl
            << "       - number of Hadamard gates     : " << nH 
            << std::endl
            << "     * Compression setting            : " << compr
            << std::endl
            << "       - CNOT compression ratio [%]   : " 
            << (1.0 - (double(nCNOT)/maxgates))*100
            << std::endl
            << "       - Ry compression ratio [%]     : " 
            << (1.0 - (double(nRY)/maxgates))*100
            << std::endl << std::endl ;

  // Transform data back to binary order, inverse FWHT, and grayscale map
  std::cout << " * Converting compressed data to grayscale image.. " 
            << std::endl ;

  qpixl::frqi::invGrayPermutation( data ) ;
  qpixl::frqi::isfwht( data ) ; 
  std::for_each( data.begin(), data.end(), []( double& n ){ n /= 2; } ) ;       
  qpixl::frqi::convertToGrayscale( data , maxval ) ;

  // Write image data to PGM
  std::cout << " * Writing compressed image to file.. " << std::endl ;
  std::string fileNamePGM ;
  fileNamePGM.append( fileNameOut ) ;
  fileNamePGM.append( ".pgm" ) ;
  s = qpixl::writePGMA(fileNamePGM , data , nrows , ncols , maxval ,
                       zeroPad) ;

  if ( s != 0 ) {
    std::cerr << "ERROR -- writing compressed image file" << std::endl;
    return -3 ;
  }

  // Simulate the circuit
  if ( sim == 1 ) { 

    // Initial state vector
    std::vector< double > state( (size_t(1)<<circuit.nbQubits()) , 0.0 );
    state[0] = 1.0 ;
    
    // Simulate
    std::cout << " * Simulating the compressed FRQI circuit.. " << std::endl ;
    circuit.simulate( state ) ;
    
    // Recover simulated image data from state vector
    for ( int i = 0; i < state.size(); i+=2 ) {
      state[i/2] = std::atan2( state[i+1] , state[i] ) ;
    }
    state.resize( state.size()/2 ) ;
    qpixl::frqi::convertToGrayscale( state , maxval ) ;

    // Write simulated image data to PGM
    std::cout << " * Writing simulated compressed image to file.. " 
              << std::endl ;
    fileNameOut.append( "_sim.pgm" ) ;
    s = qpixl::writePGMA(fileNameOut , state , nrows , ncols , maxval , 
                         zeroPad) ;

    if ( s != 0 ) {
      std::cerr << "ERROR -- writing simulated image file" << std::endl;
      return -4 ;
    }
  }

  // successful
  return 0 ;    
    
}
