#!/bin/sh -e
methods="direct"
gendir=hammersley_generator

cd hammersley_generator && ./execute_generator.sh && cd -

for testcase in $gendir/*.noref.xml.gz; do
    for method in $methods; do
	echo "======================================="
	echo "== $method"
        filename=$(basename $testcase .noref.xml.gz)

        ./scafacos_test -o $gendir/$filename.xml.gz $method $testcase
        gunzip $gendir/$filename.xml.gz

        sed 's/error_potential="1"/error_potential="1e-16"/' $gendir/$filename.xml > $gendir/tmp_file.xml
        sed 's/error_field="1"/error_field="1e-16"/' $gendir/tmp_file.xml > $gendir/$filename.xml
        rm $gendir/tmp_file.xml

        gzip $gendir/$filename.xml
	ec=$?
	if test $ec = 0; then
	    echo "PASS (exit code $ec)"
	else
	    echo "FAIL (exit code $ec)"
	fi
	echo
    done
done
