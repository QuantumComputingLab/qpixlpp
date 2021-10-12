# QPIXL++ - Quantum Image Pixel Library [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5557893.svg)](https://doi.org/10.5281/zenodo.5557893)

<p align="center"><img src="doc/doxygen/images/logo200x200.png?raw=true" /></p>

The QPIXL library currently supports the compilation of compressed quantum circuits for Flexible Representation of Quantum Images (FRQI) that contain quadratically fewer gates than previous implementations of FRQI and are NISQ-friendly as the circuits only contain CNOTs and single qubit rotation gates. 


## How to run?

The QPIXL++ package uses the CMake build system (CMake version ≥ 3.16).
The recommended way of building QPIXL++ is as follows:

1. Install

        git clone https://github.com/QuantumComputingLab/qpixlpp.git

2. CMake

        cd qpixlpp
        mkdir release
        cd release
        cmake -DCMAKE_BUILD_TYPE=Release ..
        make -j8

3. Run tests

        ./test/qpixl_tests

4. Examples

        ./examples/compressedFRQI ../examples/Example0.pgm ../examples/output 0 0

   For help

        ./examples/compressedFRQI --help

5. Generate documentation

        doxygen doxygen.dox

## References
The QPIXL++ package is based on:
- [Quantum pixel representations and compression for N-dimensional images](https://arxiv.org/abs/2110.04405),
  Mercy Amankwah, Daan Camps, E. Wes Bethel, Roel Van Beeumen, and Talita Perciano (2021)


## Developers and Contributors - Lawrence Berkeley National Laboratory
- [Daan Camps](http://campsd.github.io/) - dcamps@lbl.gov
- [Mercy Amankwah](https://mathstats.case.edu/student/mercy-amankwah/) - mercy.amankwah@case.edu<sup>1</sup>
- [E. Wes Bethel](https://dav.lbl.gov/~wes/) - ewbethel@lbl.gov
- [Talita Perciano](https://tperciano.wixsite.com/home) - tperciano@lbl.gov
- [Roel Van Beeumen](http://www.roelvanbeeumen.be/) - rvanbeeumen@lbl.gov

<sup>1</sup>Mercy Amankwah was a summer intern at Lawrence Berkeley National Laboratory during this project.

## Funding
The QPIXL++ project is supported by the Laboratory Directed Research and
Development Program of Lawrence Berkeley National Laboratory under U.S.
Department of Energy Contract No. DE-AC02-05CH11231.


## About
Quantum Image Pixel Library (QPIXL++) Copyright (c) 2021, The 
Regents of the University of California, through Lawrence Berkeley 
National Laboratory (subject to receipt of any required approvals 
from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.
       