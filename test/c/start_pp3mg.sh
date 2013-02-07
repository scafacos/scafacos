#! /bin/sh

. ../defs || exit 1

start_mpi_job -np 2 ./interface_test pp3mg ../inp_data/fortran_interface_test_inp_data.dat 1
