#! /bin/sh

# Script to generate Fortran 2003 wrappers for FFTW's MPI functions.  This
# is necessary because MPI provides no way to deal with C MPI_Comm handles
# from Fortran (where MPI_Comm == integer), but does provide a way to
# deal with Fortran MPI_Comm handles from C (via MPI_Comm_f2c).  So,
# every FFTW function that takes an MPI_Comm argument needs a wrapper
# function that takes a Fortran integer and converts it to MPI_Comm.

# pnfft.h depends on fftw3-mpi.h, fftw3.h and pfft.h
# set these paths such that the preprocessor can find the required headers
FFTW_INC=$HOME/local/pfft-1.0.6-alpha/include
PFFT_INC=$HOME/local/fftw-3.3.3/include

echo "/* Generated automatically.  DO NOT EDIT! */"
echo

echo "#include \"pnfft.h\""
echo "#include \"ipnfft.h\""
echo

# Declare prototypes using FFTW_EXTERN, important for Windows DLLs
mpicc -E pnfft.h -I${FFTW_INC} -I${PFFT_INC} |grep "pnfftl_init" |tr ';' '\n' | grep "MPI_Comm" | perl genf03-wrap.pl | grep "MPI_Fint" | sed 's/^/PNFFT_EXTERN /;s/$/;/'
mpicc -E pnfft.h -I${FFTW_INC} -I${PFFT_INC} |grep "pnfftl_init" |tr ';' '\n' | grep "MPI_Comm" | perl genf03-wrap.pl

