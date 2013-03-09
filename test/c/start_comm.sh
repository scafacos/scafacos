#! /bin/sh

. ../defs || exit 1

np=4

[ -n "$NP" ] && np=$NP

methods="direct ewald fmm memd mmm1d mmm2d p2nfft p3m pepc pp3mg vmg"

echo "--------------------------------------------------"

rm -f test_comm.output

for m in ${methods} ; do
  start_mpi_job -np $np ./test_comm $m >> test_comm.output 2>&1
  ret=$?
  echo -n "comm_test: method '$m' "
  case "$ret" in
    "0") echo "is fine!" ;;
    "1") echo "FAILED due to timeout!" ;;
    "10") echo "is not available!" ;;
    "11") echo "FAILED on 'fcs_init'" ;;
    "12") echo "FAILED on 'fcs_set_common'" ;;
    "13") echo "FAILED on 'fcs_run'" ;;
    "14") echo "FAILED on 'fcs_destroy'" ;;
    *) echo "FAILED with return code '$ret'" ;;
  esac
done

rm -f test_comm.output

echo "--------------------------------------------------"

exit 0
