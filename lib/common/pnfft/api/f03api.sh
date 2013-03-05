#! /bin/sh

# Script to generate Fortran 2003 interface declarations for FFTW from
# the fftw3.h header file.

# This is designed so that the Fortran caller can do:
#   use, intrinsic :: iso_c_binding
#   implicit none
#   include 'fftw3.f03'
# and then call the C FFTW functions directly, with type checking.

echo "! Generated automatically.  DO NOT EDIT!"
echo

# C_FFTW_R2R_KIND is determined by configure and inserted by the Makefile
# echo "  integer, parameter :: C_FFTW_R2R_KIND = @C_FFTW_R2R_KIND@"

# Extract constants
# echo "! integers"
# perl -pe 's/([A-Z0-9_]+)=([+-]?[0-9]+)/\n  integer\(C_INT\), parameter :: \1 = \2\n/g' < pnfft.h | grep 'integer(C_INT)'
# echo "! unsigned"
# perl -pe 's/#define +([A-Z0-9_]+) +\(([+-]?[0-9]+)U?\)/\n  integer\(C_INT\), parameter :: \1 = \2\n/g' < pnfft.h | grep 'integer(C_INT)'
# echo "! shifted unsigned"
# perl -pe 'if (/#define +([A-Z0-9_]+) +\(([0-9]+)U? *<< *([0-9]+)\)/) { print "\n  integer\(C_INT\), parameter :: $1 = ",$2 << $3,"\n"; }' < pnfft.h | grep 'integer(C_INT)'
# echo "! redirections"
# perl -pe 'if (/#define +([A-Z0-9_]+) +\(\(([A-Z0-9_| ]+)\)\)/) { print "\n  integer\(C_INT\), parameter :: $1 = $2\n"; }' < pnfft.h | grep 'integer(C_INT)' | sed 's/| / + /g'

# Extract function declarations
for p in $*; do
    if test "$p" = "d"; then p=""; fi

    echo
    cat <<EOF
EOF

    echo
    echo "  interface"
#     gcc -E pnfft.h |grep "pnfft${p}_create_procmesh" |tr ';' '\n'
    gcc -E pnfft.h -I$HOME/local/pfft-1.0.6-alpha/include -I$HOME/local/fftw-3.3.3/include |grep "pnfft${p}_init" |tr ';' '\n'

#     gcc -E pnfft.h |grep "pnfft${p}_create_procmesh" |tr ';' '\n' | perl genf03.pl
#     gcc -E pnfft.h |grep "pnfft${p}_create_procmesh" |tr ';' '\n' | grep -v "pnfft${p}_trafo(" | perl genf03.pl
    echo "  end interface"

done
