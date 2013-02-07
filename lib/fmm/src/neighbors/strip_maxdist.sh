#!/bin/bash

echo "Stripping constants to file: include/fmmdist_maxdist.h"
#grepping, cutting, space removal, add comma, remove last line comma
grep " maxdist(" get_maxdist.compat.f90 | colrm 1 25 | sed 's/^\s*//;s/$/,/;$s/,//' > include/fmmdist_maxdist.h

