#! /bin/sh

. ../defs || exit 1

start_mpi_job -np 1 ./test_p2nfft_nonperiodic 
