#!/bin/bash

echo "Extract data from getdist.f90"
grep distx\(..[0-9],..[0-9]\)= getdist.f90 | colrm 1 24 | sed '''s/ //' | sed '''s/$/,/' > include/distx.h
grep disty\(..[0-9],..[0-9]\)= getdist.f90 | colrm 1 24 | sed '''s/ //' | sed '''s/$/,/' > include/disty.h
grep distz\(..[0-9],..[0-9]\)= getdist.f90 | colrm 1 24 | sed '''s/ //' | sed '''s/$/,/' > include/distz.h

echo "Extract data from getdistms.f90"
grep maxdist\(..[0-9],..[0-9]\)= getdistms.f90 | colrm 1 27 | sed '''s/ //' | sed '''s/$/,/' > include/maxdist.h
grep dist3dsquared\(..[0-9],..[0-9]\)= getdistms.f90 | colrm 1 35 | sed '''s/ //' | sed '''s/$/,/' > include/dist3dsquared.h
