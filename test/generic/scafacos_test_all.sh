#!/bin/sh
methods="direct ewald fmm memd mmm1d mmm2d p2nfft p3m pepc pp3mg vmg"
testcases="systems/3d-periodic/cloud-wall.xml.gz"

for testcase in $testcases; do
    for method in $methods; do
	echo "======================================="
	echo "== $method"
	./scafacos_test $method $testcase
	ec=$?
	if test $ec = 0; then
	    echo "PASS (exit code $ec)"
	else
	    echo "FAIL (exit code $ec)"
	fi
	echo
    done
done
