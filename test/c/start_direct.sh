#! /bin/sh

#IGNORE_RUNCHECKS=yes

. ../defs || exit 1

np=2

[ -n "$NP" ] && np=$NP

start_mpi_job -np $np ./test_direct ../inp_data/direct/direct_conf.in
