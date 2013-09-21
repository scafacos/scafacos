#!/bin/sh
#set -xv

SL_CONFIG=../../common/sl/support/sl_config.sh

${SL_CONFIG} --dst-sl=. --source-ref=../../common/sl --am-libname=libsl_pepc.a --extra-prefix=fcs_pepc_ --config-not=pepcparts
