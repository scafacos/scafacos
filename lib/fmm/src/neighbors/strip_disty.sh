#!/bin/bash

echo "Stripping constants to file: include/fmmdist_disty.h"
#grepping, cutting, space removal, add comma, remove last line comma
grep " disty(" get_disty.compat.f90 | colrm 1 23 | sed 's/^\s*//;s/$/,/;$s/,//' > include/fmmdist_disty.h
