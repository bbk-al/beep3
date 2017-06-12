#!/bin/bash
# Build old version using old compiler and libgsl
# Don't assume profile has not been run!!
export LD_LIBRARY_PATH=
export PATH=${PATH#/d/mw6/u/la002/pjt/build/bin:}
g++ -I../../../../build/include tsphmb.cpp -L../../../../build/gsl/gsl-1.12/.libs -lgsl -L../../../../build/gsl/gsl-1.12/cblas/.libs -lgslcblas -o tsphmb

# Build new version using new compiler and libgsl
. ../../../../profile
g++ -I../../../../build/include tsphmt.cpp -L../../../../build/gsl/gsl-2.3/.libs -lgsl -L../../../../build/gsl/gsl-2.3/cblas/.libs -lgslcblas -o tsphmt

# Run old and new versions

LD_LIBRARY_PATH=../../../../build/gsl/gsl-1.12/.libs:../../../../build/gsl/gsl-1.12/cblas/.libs ./tsphmb | sort -k1,2 >tsphmb.out 2>&1
./tsphmt | sort -k1,2 >tsphmt.out 2>&1

diff tsphm?.out
