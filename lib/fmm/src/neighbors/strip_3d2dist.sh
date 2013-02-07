#!/bin/bash

echo "Stripping constants to file: include/fmmdist_3d2dist.h"
#grepping, cutting, space removal, add comma, remove last line comma
grep " 3d2dist(" get_3d2dist.compat.f90 | colrm 1 25 | sed 's/^\s*//;s/$/,/;$s/,//' > include/fmmdist_3d2dist.h

