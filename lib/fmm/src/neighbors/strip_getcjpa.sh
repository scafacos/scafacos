#!/bin/bash

echo "Extract data from getcjpa.f90"
grep casejump\(..[0-9]\) getcjpa.f90 | colrm 1 23 | sed '''s/ //' | sed '''s/$/,/' > include/casejump.h
