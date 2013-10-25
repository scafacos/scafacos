#!/bin/sh
#set -xv

SL_CONFIG=../../common/sl/support/sl_config.sh

${SL_CONFIG} --dst-sl=. --source-ref=../../common/sl --am-libname=libsl_pepc.a --am-libprefix=fcs_pepc_ --extra-prefix=fcs_pepc_ --config-not=pepcparts
