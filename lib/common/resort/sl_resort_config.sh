#!/bin/sh
#set -xv

SL=../sl

SL_CONFIG=${SL}/support/sl_config.sh

${SL_CONFIG} --src-sl="${SL}" --dst-sl=. --source=single --source-link --am-libname=libfcs_resort.a --am-libprefix=fcs_resort_ --extra-prefix=fcs_resort_

sed -i -e 's!^noinst_LTLIBRARIES\(.*libfcs_resort.la\)$!if ENABLE_SINGLE_LIB\nnoinst_LTLIBRARIES\1\nelse\nlib_LTLIBRARIES\1\n\endif!' \
       -e 's!\(-I\$(srcdir)/include\)!-I$(top_srcdir)/lib \1!' \
       Makefile.am
