#! /bin/sh

. ../defs || exit 1

start_mpi_job -np 1 ./test_vmg_periodic
#start_mpi_job -np 1 ./test_vmg_open
