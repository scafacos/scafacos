#!/bin/bash

echo "Stripping constants to file: include/fmmdist_distz.h"
#grepping, cutting, space removal, add comma, remove last line comma
grep " distz(" get_distz.compat.f90 | colrm 1 23 | sed 's/^\s*//;s/$/,/;$s/,//' > include/fmmdist_distz.h
