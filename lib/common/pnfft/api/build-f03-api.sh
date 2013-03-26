#! /bin/sh
## Script to generate Fortran 2003 interface and wrappers. 

sh f03-api.sh d f > pnfft.f03.in
sh f03-api.sh l > pnfftl.f03.in
sh f03-wrap.sh > f03-wrap.c

