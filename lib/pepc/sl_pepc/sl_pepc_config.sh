#!/bin/sh
#set -xv

SL=../../common/sl

SL_CONFIG=${SL}/support/sl_config.sh

${SL_CONFIG} --src-sl="${SL}" --dst-sl=. --source=single --source-link --am-libname=libsl_pepc.a --am-libprefix=fcs_pepc_ --extra-prefix=fcs_pepc_ --config-not=pepcparts
