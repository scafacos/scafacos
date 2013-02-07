#!/bin/bash

echo "Stripping constants to file: include/fmmdist_distx.h"
#grepping, cutting, space removal, add comma, remove last line comma
grep " distx(" get_distx.compat.f90 | colrm 1 23 | sed 's/^\s*//;s/$/,/;$s/,//' > include/fmmdist_distx.h
