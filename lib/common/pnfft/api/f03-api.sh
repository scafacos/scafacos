#! /bin/sh

# Script to generate Fortran 2003 interface declarations for FFTW from
# the fftw3.h header file.

# This is designed so that the Fortran caller can do:
#   use, intrinsic :: iso_c_binding
#   implicit none
#   include 'fftw3.f03'
# and then call the C FFTW functions directly, with type checking.


# pnfft.h depends on fftw3-mpi.h, fftw3.h and pfft.h
# set these paths such that the preprocessor can find the required headers
FFTW_INC=$HOME/local/pfft-1.0.6-alpha/include
PFFT_INC=$HOME/local/fftw-3.3.3/include

echo "! Generated automatically.  DO NOT EDIT!"

# C_FFTW_R2R_KIND is determined by configure and inserted by the Makefile
# echo "  integer, parameter :: C_FFTW_R2R_KIND = @C_FFTW_R2R_KIND@"

# Extract constants
echo
echo "! integers"
perl -pe 's/([A-Z0-9_]+)=([+-]?[0-9]+)/\n  integer\(C_INT\), parameter :: \1 = \2\n/g' < pnfft.h | grep 'integer(C_INT)'
echo
echo "! unsigned"
perl -pe 's/#define +([A-Z0-9_]+) +\(([+-]?[0-9]+)U?\)/\n  integer\(C_INT\), parameter :: \1 = \2\n/g' < pnfft.h | grep 'integer(C_INT)'
echo
echo "! shifted unsigned"
perl -pe 'if (/#define +([A-Z0-9_]+) +\(([0-9]+)U? *<< *([0-9]+)\)/) { print "\n  integer\(C_INT\), parameter :: $1 = ",$2 << $3,"\n"; }' < pnfft.h | grep 'integer(C_INT)'
echo
echo "! redirections"
perl -pe 'if (/#define +([A-Z0-9_]+) +\(\(([A-Z0-9_| ]+)\)\)/) { print "\n  integer\(C_INT\), parameter :: $1 = $2\n"; }' < pnfft.h | grep 'integer(C_INT)' | sed 's/| / + /g'

# Extract function declarations
for p in $*; do
    if test "$p" = "d"; then p=""; fi

    echo
#     cat <<EOF
#   type, bind(C) :: fftw${p}_mpi_ddim
#      integer(C_INTPTR_T) n, ib, ob
#   end type fftw${p}_mpi_ddim   cat <<EOF
# EOF

    echo
    echo "  interface"
    grep -v 'mpi.h' pnfft.h | gcc -E - -I${FFTW_INC} -I${PFFT_INC} |grep "pnfft${p}_init" |tr ';' '\n' | perl genf03-api.pl
    echo "  end interface"

done
