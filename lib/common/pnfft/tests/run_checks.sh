#!/bin/sh -e


## Extra checks:
test_files_4="$test_files simple_test simple_test_c2r"
test_files_4="$test_files check_trafo_2d check_trafo_transposed_2d"

## 3D:
test_files_2="$test_files check_trafo_2d"

## 3D on 3D procmesh:
test_files_8="$test_files pnfft_test pnfft_test_adv"
test_files_8="$test_files check_trafo check_trafo_grad"
test_files_8="$test_files check_trafo_transposed check_trafo_transposed_grad"
test_files_8="$test_files check_trafo_vs_ndft check_trafo_vs_ndft_c2r check_trafo_grad_vs_ndft check_trafo_vs_ndft_transposed_2d"

test_files_8="$test_files check_adj check_adj_transposed"
test_files_8="$test_files check_adjvs_ndft check_adj_vs_ndft_c2r"
test_files_8="$test_files pnfft_test pnfft_test_adv"

test_files_8="$test_files simple_test_c2r_c2c_compare_complex.c"
test_files_8="$test_files simple_test_c2r_c2c_compare_grad.c"
test_files_8="$test_files simple_test_c2r_c2c_compare_real.c"
test_files_8="$test_files simple_test_c2r_c2c_compare_timer.c"

testpath=`dirname ${0}`
 
## Run tests with 2 processes
for name in $test_files_2; do
  echo "## mpirun -np 2 ${name}:"
  mpirun -np 2 ${testpath}/${name}
done

## Run tests with 4 processes
for name in $test_files_4; do
  echo "## mpirun -np 4 ${name}:"
  mpirun -np 4 ${testpath}/${name}
done

## Run tests with 8 processes
for name in $test_files_8; do
  echo "## mpirun -np 8 ${name}:"
  mpirun -np 8 ${testpath}/${name}
done
