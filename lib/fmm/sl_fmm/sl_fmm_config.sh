#!/bin/sh
#set -xv

SL=../../common/sl

SL_CONFIG=${SL}/support/sl_config.sh

${SL_CONFIG} --src-sl="${SL}" --dst-sl=. --source=single --source-link --am-libname=libsl_fmm.a --am-libprefix=fcs_fmm_ --extra-prefix=fcs_fmm_
