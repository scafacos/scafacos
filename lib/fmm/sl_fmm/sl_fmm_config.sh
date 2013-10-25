#!/bin/sh
#set -xv

SL_CONFIG=../../common/sl/support/sl_config.sh

${SL_CONFIG} --dst-sl=. --source-ref=../../common/sl --am-libname=libsl_fmm.a --am-libprefix=fcs_fmm_ --extra-prefix=fcs_fmm_
