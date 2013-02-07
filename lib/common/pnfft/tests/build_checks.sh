#!/bin/sh -e

FCSDIR=$HOME/Repos/scafacos
FCSCOMMONDIR=$FCSDIR/trunk/lib/common

FFTWDIR=$HOME/local/fftw-3.3-debug
PFFTLIBDIR=$FCSCOMMONDIR/pfft/tmp_build_debug/.libs
PFFTINCDIR=$FCSCOMMONDIR/pfft/api


SRCDIR=../..

## Extra checks:
test_files="$test_files simple_test"
test_files="$test_files pnfft_check"
test_files="$test_files pnfft_check_2d"
test_files="$test_files pnfft_check_grad"

# Cecam workshop
test_files="$test_files pnfft_test"
test_files="$test_files pnfft_test_adv"

PFFTLIB=$PFFTLIBDIR/lib*pfft*.a
PNFFTLIB=../.libs/lib*pnfft*.a

cd .. && make && cd - 
for name in $test_files; do
  mpicc $SRCDIR/tests/${name}.c -o ${name} \
    -I$FFTWDIR/include -I$SRCDIR/api -I$PFFTINCDIR -std=gnu99 -O0 -g -Wall \
    $PNFFTLIB $PFFTLIB $FFTWDIR/lib/libfftw3_mpi.a $FFTWDIR/lib/libfftw3.a
done
