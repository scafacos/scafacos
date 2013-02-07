#! /bin/sh

. ../defs || exit 1

start_mpi_job -np 2 ./test_pepc
