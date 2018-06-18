#!/bin/sh
#set -xv

SL=../sl

SL_CONFIG=${SL}/support/sl_config.sh

${SL_CONFIG} --src-sl="${SL}" --dst-sl=. --source=single --source-link --am-libname=libfcs_near.a --am-libprefix=fcs_near_ --extra-prefix=fcs_near_

sed -i -e 's!^noinst_LTLIBRARIES\(.*libfcs_near.la\)$!if ENABLE_SINGLE_LIB\nnoinst_LTLIBRARIES\1\nelse\nlib_LTLIBRARIES\1\n\endif!' \
       -e 's!\(-I\$(srcdir)/include\)!-I$(top_srcdir)/lib \1!' \
       -e '/sl_sub_libs =/{s! \([^ ]*_dip_.*_dip_[^ ]*\)\(.*\)$!\2\nif ENABLE_DIPOLES\nsl_sub_libs += \1\nendif\n!}' \
       Makefile.am
