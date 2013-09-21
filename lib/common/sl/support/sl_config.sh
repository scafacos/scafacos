#!/bin/sh

#  
#  Copyright (C) 2011, 2012, 2013 Michael Hofmann
#  
#  This file is part of ScaFaCoS.
#  
#  ScaFaCoS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ScaFaCoS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  

#  
#  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
#  

#set -xv

VERSION=0.1

# make an own PWD that is always ending with /
my_pwd=${PWD%/}/

# script name
me=${0##*/}

# script directory
me_dir=${0%/*}
[ "${me_dir}" != "${0}" ] && me_dir="${me_dir}/"

# default config file
cfg_file="${me_dir}sl_config.cfg"
[ -f "${cfg_file}" ] && . "${cfg_file}"


# set some defaults
dstdir_prefix=sl_
header_prefix=sl_

src_sl=${me_dir#./}
src_sl_src=src/
src_sl_src_sub=true
src_sl_extra=extra/
src_sl_extra_sub=true
src_sl_scripts=scripts/
src_sl_scripts_sub=true
src_sl_adds=src/
src_sl_adds_sub=true

dst_sl=dst/
dst_sl_src=src/
dst_sl_src_sub=true
dst_sl_scripts=scripts/
dst_sl_scripts_sub=true
dst_sl_extra=extra/
dst_sl_extra_sub=true
dst_sl_adds=src/
dst_sl_adds_sub=true

dst_if=include/
dst_if_sub=true

cfg_config=
cfg_config_not=

cfg_quiet=
cfg_verbose=
cfg_debug=

cfg_local=

cfg_far=perl
cfg_no_files=
cfg_testing=
cfg_di_versions=false  # false|config|makefile
cfg_interface=true
cfg_makefile=true
cfg_makefile_in=
cfg_prefix=
cfg_autoconf=
cfg_automake=
cfg_automake_libname=libsl.a
cfg_automake_noinst=
cfg_source=
cfg_source_rename=
cfg_source_ref=
cfg_source_ref_set=
cfg_extra=true
cfg_extra_prefix=
cfg_extra_have_h="zmpi_local zmpi_tools zmpi_ataip zmpi_atasp"
cfg_extra_include="z_pack zmpi_local zmpi_tools"
cfg_extra_wrap_have_h="zmpi_local zmpi_tools zmpi_ataip zmpi_atasp"
cfg_extras="z_pack zmpi_local zmpi_tools zmpi_ataip zmpi_atasp"
cfg_scripts=
cfg_tuning=

cfg_license_file=

cfg_makefile_makefile_in=
cfg_makefile_target=
cfg_makefile_sl_use_mpi=
cfg_makefile_all=libsl
cfg_makefile_quiet=
cfg_makefile_fixed=
cfg_makefile_bulk_ar=

cfg_makefile_incdir=

cfg_makefile_libdir=

cfg_makefile_wrapper_src=
cfg_makefile_wrapper_prefix=
cfg_makefile_wrapper_mpi_src=
cfg_makefile_wrapper_mpi_prefix=mpi_

cfg_makefile_exec="a.out"
cfg_makefile_exec_src=
cfg_makefile_exec_prefix=
cfg_makefile_exec_mpi_src=
cfg_makefile_exec_mpi_prefix=mpi_

cmd_sed="sed"
cmd_cat="cat"
cmd_grep="grep"
cmd_cp="cp"
#cmd_cpp="cpp"
cmd_mkdir="mkdir -p"

[ -f "${me_dir}sl_config.cfg" ] && . "${me_dir}sl_config.cfg"

# default file descriptors
exec 3>&1         # 3 = progress, default to stdout
exec 4>/dev/null  # 4 = verbose, default to /dev/null
exec 5>/dev/null  # 5 = debug, default to /dev/null


pout()
{
  local n="\n"

  [ "$1" = "-n" ] && { n="" ; shift ; }

  printf "$*${n}" >&1
}

perror()
{
  local n="\n"

  [ "$1" = "-n" ] && { n="" ; shift ; }

  printf "${me}: ERROR: $*${n}" >&2
}

perror_exit()
{
  local x=1

  [ "$1" = "-x" ] && { x="$2" ; shift 2 ; }

  perror $*

  exit ${x}
}

pprogress()
{
  local n="\n"

  [ "$1" = "-n" ] && { n="" ; shift ; }

  printf "$*${n}" >&3
}

pverbose()
{
  local n="\n"

  [ "$1" = "-n" ] && { n="" ; shift ; }

  printf "$*${n}" >&4
}

pdebug()
{
  local n="\n"

  [ "$1" = "-n" ] && { n="" ; shift ; }

  printf "${me}: DEBUG: $*${n}" >&5
}

print_settings()
{
  local cmd_print="$1"
  local prefix="$2"

  local list1="my_pwd cfg_config cfg_config_not"
  local list2="cfg_di_versions cfg_makefile_sl_use_mpi cfg_local cfg_makefile cfg_makefile_in cfg_prefix cfg_autoconf cfg_automake cfg_interface cfg_source cfg_source_rename cfg_source_ref cfg_source_ref_set"

  local list_src="src_sl src_sl_src src_sl_src_sub src_sl_extra src_sl_extra_sub src_sl_scripts src_sl_scripts_sub src_sl_adds src_sl_adds_sub"
  local list_dst="dst_sl dst_sl_src dst_sl_src_sub dst_sl_extra dst_sl_extra_sub dst_sl_scripts dst_sl_scripts_sub dst_sl_adds dst_sl_adds_sub dst_if dst_if_sub"
  local list_dst="ref_sl ref_sl_src ref_sl_extra ref_sl_scripts ref_sl_adds"

  local list_mf="cfg_makefile_makefile_in cfg_makefile_target cfg_makefile_sl_use_mpi cfg_makefile_all cfg_makefile_quiet cfg_makefile_fixed cfg_makefile_bulk_ar cfg_makefile_incdir cfg_makefile_libdir"
  local list_mf_wrap="cfg_makefile_wrapper_src cfg_makefile_wrapper_prefix cfg_makefile_wrapper_mpi_src cfg_makefile_wrapper_mpi_prefix"
  local list_mf_exec="cfg_makefile_exec cfg_makefile_exec_src cfg_makefile_exec_prefix cfg_makefile_exec_mpi_src cfg_makefile_exec_mpi_prefix"

  ${cmd_print} ""
  ${cmd_print} "${prefix}settings:"

  for s in ${list1} ${list2} ${list_src} ${list_dst} ${list_mf} ${list_mf_wrap} ${list_mf_exec} ; do
    eval v=\$${s}
    ${cmd_print} "${prefix} ${s}: ${v}"
  done
}

print_commands()
{
  local cmd_print="$1"
  local prefix="$2"

  local list1="cmd_cat cmd_grep cmd_sed cmd_cp cmd_mkdir cmd_far"

  ${cmd_print} ""
  ${cmd_print} "${prefix}commands:"

  for s in ${list1} ; do
    eval v=\$${s}
    ${cmd_print} "${prefix} ${s}: ${v}"
  done
}

canonize_dir()
{
  local xdir
  local v
  local vnew

  while [ -n "$1" ] ; do
    xdir=$1
    shift
    eval v=\$${xdir}
#    pdebug "canonize ${xdir}: $v"
    vnew=
    for d in ${v} ; do
      eval d=`echo ${d}`
      [ "${d%/}" = "${d}" ] && d="${d}/"
      vnew="${vnew} ${d}"
    done
    eval ${xdir}=\"${vnew# }\"
#    pdebug "canonized ${xdir}: ${vnew}"
  done
}

mkdir_save()
{
  while [ -n "$*" ] ; do
    [ -n "$1" -a ! -d "$1" ] && ${cmd_mkdir} $1
    shift
  done
}


shift_opt=

#check_arg_dir()
#{
#  local key=$1
#  local var=$2
#  local opt=$3
#  local val=$4
#
#  echo "'$key' / '$var' / '$opt' / '$arg'"
#
#  if [ "-${var}" = "${opt}" -o "--${var}" = "${opt}" ] ; then
#
#    pdebug "option: $key = ${optarg}"
#
#    eval $key="${optarg}"
#    
#    shift_opt="1"
#
#  else
#
#    shift_opt=
#
#  fi
#}
#
#all=
#
#while [ -n "$1" ] ; do
#
#  opt="$1"
#  optarg=`expr "x$1" : 'x[^=]*=\(.*\)'`
#  opt="${opt%%=*}"
#  
#  if [ "${opt}" = "$1" ] ; then
#    shift
#    optarg=$1
#  fi
#
#  echo "opt: '${opt}' / '${optarg}'"
#
#  check_arg_dir "src_sl" "src-sl" ${opt} ${optargs}
#
#  echo "shift_opt: ${shift_opt}"
#
#  [ -n "${shift_opt}" ] && shift ${shift_opt}
#
#  all="${all} ${opt}"
#
#done



# process all arguments
prev=
append=
all=$*
for opt in $all ; do

  if [ -n "${prev}" ] ; then

    opt_isset=true

  else

    opt_isset=

    optarg=`expr "x${opt}" : 'x[^=]*=\(.*\)'`

  case ${opt} in
    -h | --help)
      pout " help is coming soon!"
      exit 0
      ;;
    -V | --version)
      pout " this is version v${VERSION}"
      exit 0
      ;;
    -q | --quiet | --silent)
      cfg_quiet=true
      exec 3> /dev/null
      ;;
    -v | -verbose | --verbose | --verbos | --verbo | --verb)
      cfg_verbose=true
      exec 4>&1
      ;;
    -d | --debug)
      cfg_debug=true
      exec 5>&1
      ;;

    # find-and-replace: sed, perl (default)
    -far | --far)
      prev=cfg_far
      ;;
    -far=* | --far=*)
      cfg_far=${optarg}
      ;;

    -no-files | --no-files)
      cfg_no_files=true  # FIXME: not so easy, what should I do with 'echo ... > ...'?
      ;;
    -testing | --testing)
      cfg_testing=true
      ;;

    -config-devel | --config-devel)
      dst_sl=
      cfg_local=true
      cfg_source=ref
      cfg_makefile=
      cfg_autoconf=
      cfg_interface=
      cfg_scripts=
      cfg_extra=true
      cfg_extra_prefix=
      ;;

    -local | --local)
      cfg_local=true
      ;;

    -clean | --clean)
      cfg_clean=true
      ;;


    -di-versions | --di-versions)
     prev=cfg_di_versions
     ;;
    -di-versions=* | --di-versions=*)
     cfg_di_versions=${optarg}
     ;;
    -interface | --interface)
      cfg_interface=true
      ;;
    -no-interface | --no-interface)
      cfg_interface=
      ;;
    -makefile | --makefile)
      cfg_makefile=true
      ;;
    -no-makefile | --no-makefile)
      cfg_makefile=
      ;;
    -makefile-in | --makefile-in)
      prev=cfg_makefile_in
      append=true
      ;;
    -makefile-in=* | --makefile-in=*)
      prev=cfg_makefile_in
      append=true
      opt=${optarg}
      opt_isset=true
      ;;
    -autoconf | --autoconf)
      cfg_autoconf=true
      ;;
    -no-autoconf | --no-autoconf)
      cfg_autoconf=
      ;;
    -automake | --automake)
      cfg_automake=true
      ;;
    -no-automake | --no-automake)
      cfg_automake=
      ;;
    -source | --source | -src | --src) # ref | single | sep[arate] | <empty> = ref
      prev=cfg_source
      ;;
    -source=* | --source=* | -src=* | --src=*) # ref | single | sep[arate] | <empty> = ref
      prev=cfg_source
      opt=${optarg}
      opt_isset=true
      ;;
    -source-rename | --source-rename | -src-rename | --src-rename)
      cfg_source_rename=true
      ;;
    -source-ref | --source-ref | -src-ref | --src-ref)
      prev=cfg_source_ref
      cfg_source_ref_set=true
      ;;
    -source-ref=* | --source-ref=* | -src-ref=* | --src-ref=*)
      prev=cfg_source_ref
      opt=${optarg}
      opt_isset=true
      cfg_source_ref_set=true
      ;;
    -prefix | --prefix)
      prev=cfg_prefix
      append=true
      ;;
    -prefix=* | --prefix=*)
      prev=cfg_prefix
      append=true
      opt=${optarg}
      opt_isset=true
      ;;
    -extra | --extra)
      cfg_extra=true
      ;;
    -no-extra | --no-extra)
      cfg_extra=false
      ;;
    -extra-prefix | --extra-prefix)
      prev=cfg_extra_prefix
      ;;
    -extra-prefix=* | --extra-prefix=*)
      cfg_extra_prefix=${optarg}
      ;;
    -scripts | --scripts)
      cfg_scripts=true
      ;;
    -no-scripts | --no-scripts)
      cfg_scripts=
      ;;
    -tuning | --tuning)
      cfg_tuning=true
      ;;
    -no-tuning | --no-tuning)
      cfg_tuning=
      ;;
    -tune-auto-dummy | --tune-auto-dummy)
      cfg_tune_auto_dummy=true
      ;;

    -license-file | --license-file)
      prev=cfg_license_file
      ;;
    -license-file=* | --license-file=*)
      cfg_license_file=${optarg}
      ;;

    -src-sl | --src-sl)
      prev=src_sl
      ;;
    -src-sl=* | --src-sl=*)
      src_sl=${optarg}
      ;;
    -src-sl-src | --src-sl-src)
      prev=src_sl_src
      src_sl_src_sub=
      ;;
    -src-sl-src=* | --src-sl-src=*)
      src_sl_src=${optarg}
      src_sl_src_sub=
      ;;
    -src-sl-src-sub | --src-sl-src-sub)
      prev=src_sl_src
      src_sl_src_sub=true
      ;;
    -src-sl-src-sub=* | --src-sl-src-sub=*)
      src_sl_src=${optarg}
      src_sl_src_sub=true
      ;;
    -src-sl-scripts | --src-sl-scripts)
      prev=src_sl_scripts
      src_sl_scripts_sub=
      ;;
    -src-sl-scripts=* | --src-sl-scripts=*)
      src_sl_scripts=${optarg}
      src_sl_scripts_sub=
      ;;
    -src-sl-scripts-sub | --src-sl-scripts-sub)
      prev=src_sl_scripts
      src_sl_scripts_sub=true
      ;;
    -src-sl-scripts-sub=* | --src-sl-scripts-sub=*)
      src_sl_scripts=${optarg}
      src_sl_scripts_sub=true
      ;;
    -src-sl-adds | --src-sl-adds)
      prev=src_sl_adds
      src_sl_adds_sub=
      ;;
    -src-sl-adds=* | --src-sl-adds=*)
      src_sl_adds=${optarg}
      src_sl_adds_sub=
      ;;
    -src-sl-adds-sub | --src-sl-adds-sub)
      prev=src_sl_adds
      src_sl_adds_sub=true
      ;;
    -src-sl-adds-sub=* | --src-sl-adds-sub=*)
      src_sl_adds=${optarg}
      src_sl_adds_sub=true
      ;;

    -config | --config)
      prev=cfg_config
      ;;
    -config=* | --config=*)
      cfg_config=${optarg}
      ;;
    -config-not | --config-not)
      prev=cfg_config_not
      ;;
    -config-not=* | --config-not=*)
      cfg_config_not=${optarg}
      ;;

    -dst-sl | --dst-sl)
      prev=dst_sl
      ;;
    -dst-sl=* | --dst-sl=*)
      dst_sl=${optarg}
      ;;
    -dst-sl-src | --dst-sl-src)
      prev=dst_sl_src
      dst_sl_src_sub=
      ;;
    -dst-sl-src=* | --dst-sl-src=*)
      dst_sl_src=${optarg}
      dst_sl_src_sub=
      ;;
    -dst-sl-src-sub | --dst-sl-src-sub)
      prev=dst_sl_src
      dst_sl_src_sub=true
      ;;
    -dst-sl-src-sub=* | --dst-sl-src-sub=*)
      dst_sl_src=${optarg}
      dst_sl_src_sub=true
      ;;
    -dst-sl-scripts | --dst-sl-scripts)
      prev=dst_sl_scripts
      dst_sl_scripts_sub=
      ;;
    -dst-sl-scripts=* | --dst-sl-scripts=*)
      dst_sl_scripts=${optarg}
      dst_sl_scripts_sub=
      ;;
    -dst-sl-scripts-sub | --dst-sl-scripts-sub)
      prev=dst_sl_scripts
      dst_sl_scripts_sub=true
      ;;
    -dst-sl-scripts-sub=* | --dst-sl-scripts-sub=*)
      dst_sl_scripts=${optarg}
      dst_sl_scripts_sub=true
      ;;
    -dst-sl-adds | --dst-sl-adds)
      prev=dst_sl_adds
      dst_sl_adds_sub=
      ;;
    -dst-sl-adds=* | --dst-sl-adds=*)
      dst_sl_adds=${optarg}
      dst_sl_adds_sub=
      ;;
    -dst-sl-adds-sub | --dst-sl-adds-sub)
      prev=dst_sl_adds
      dst_sl_adds_sub=true
      ;;
    -dst-sl-adds-sub=* | --dst-sl-adds-sub=*)
      dst_sl_adds=${optarg}
      dst_sl_adds_sub=true
      ;;

    -mf-target | --mf-target | -target | --target)
      prev=cfg_makefile_target
      ;;
    -mf-target=* | --mf-target=* | -target=* | --target=*)
      cfg_makefile_target=${optarg}
      ;;
    -mf-makefile-in | --mf-makefile-in | -makefile-in | --makefile-in)
      prev=cfg_makefile_makefile_in
      ;;
    -mf-makefile-in=* | --mf-makefile-in=* | -makefile-in=* | --makefile-in=*)
      cfg_makefile_makefile_in=${optarg}
      ;;
    -mf-sl-use-mpi | --mf-sl-use-mpi | -sl-use-mpi | --sl-use-mpi)
      prev=cfg_makefile_sl_use_mpi
      ;;
    -mf-sl-use-mpi=* | --mf-sl-use-mpi=* | -sl-use-mpi=* | --sl-use-mpi=*)
      cfg_makefile_sl_use_mpi=${optarg}
      ;;
    -mf-all | --mf-all | -all | --all)
      prev=cfg_makefile_all
      ;;
    -mf-all=* | --mf-all=* | -all=* | --all=*)
      cfg_makefile_all=${optarg}
      ;;
    -mf-quiet | --mf-quiet | -quiet-makefile | --quiet-makefile) # deprecated: quiet-makefile
      cfg_makefile_quiet=true
      ;;
    -no-mf-quiet | --no-mf-quiet)
      cfg_makefile_quiet=
      ;;
    -mf-fixed | --mf-fixed | -fixed | --fixed) # deprecated: fixed
      cfg_makefile_fixed=true
      ;;
    -no-mf-fixed | --no-mf-fixed)
      cfg_makefile_fixed=
      ;;
    -mf-bulk-ar | --mf-bulk-ar)
      cfg_makefile_bulk_ar=true
      ;;
    -no-mf-bulk-ar | --no-mf-bulk-ar)
      cfg_makefile_bulk_ar=
      ;;

    -mf-incdir | --mf-incdir)
      prev=cfg_makefile_incdir
      append=true
      ;;
    -mf-incdir=* | --mf-incdir=*)
      prev=cfg_makefile_incdir
      append=true
      opt=${optarg}
      opt_isset=true
      ;;

    -mf-libdir | --mf-libdir)
      prev=cfg_makefile_libdir
      ;;
    -mf-libdir=* | --mf-libdir=*)
      cfg_makefile_libdir=${optarg}
      ;;

    -mf-wrapper-src | --mf-wrapper-src)
      prev=cfg_makefile_wrapper_src
      ;;
    -mf-wrapper-src=* | --mf-wrapper-src=*)
      cfg_makefile_wrapper_src=${optarg}
      ;;
    -mf-wrapper-prefix | --mf-wrapper-prefix)
      prev=cfg_makefile_wrapper_prefix
      ;;
    -mf-wrapper-prefix=* | --mf-wrapper-prefix=*)
      cfg_makefile_wrapper_prefix=${optarg}
      ;;
    -mf-wrapper-mpi-src | --mf-wrapper-mpi-src)
      prev=cfg_makefile_wrapper_mpi_src
      ;;
    -mf-wrapper-mpi-src=* | --mf-wrapper-mpi-src=*)
      cfg_makefile_wrapper_mpi_src=${optarg}
      ;;
    -mf-wrapper-mpi-prefix | --mf-wrapper-mpi-prefix)
      prev=cfg_makefile_wrapper_mpi_prefix
      ;;
    -mf-wrapper-mpi-prefix=* | --mf-wrapper-mpi-prefix=*)
      cfg_makefile_wrapper_mpi_prefix=${optarg}
      ;;

    -mf-exec | --mf-exec | -exec | --exec)
     prev=cfg_makefile_exec
     ;;
    -mf-exec=* | --mf-exec=* | -exec=* | --exec=*)
     cfg_makefile_exec=${optarg}
     ;;

    -mf-exec-src | --mf-exec-src)
      prev=cfg_makefile_exec_src
      ;;
    -mf-exec-src=* | --mf-exec-src=*)
      cfg_makefile_exec_src=${optarg}
      ;;
    -mf-exec-prefix | --mf-exec-prefix)
      prev=cfg_makefile_exec_prefix
      ;;
    -mf-exec-prefix=* | --mf-exec-prefix=*)
      cfg_makefile_exec_prefix=${optarg}
      ;;
    -mf-exec-mpi-src | --mf-exec-mpi-src)
      prev=cfg_makefile_exec_mpi_src
      ;;
    -mf-exec-mpi-src=* | --mf-exec-mpi-src=*)
      cfg_makefile_exec_mpi_src=${optarg}
      ;;
    -mf-exec-mpi-prefix | --mf-exec-mpi-prefix)
      prev=cfg_makefile_exec_mpi_prefix
      ;;
    -mf-exec-mpi-prefix=* | --mf-exec-mpi-prefix=*)
      cfg_makefile_exec_mpi_prefix=${optarg}
      ;;

    -am-libname=* | --am-libname=*)
      cfg_automake_libname=${optarg}
      ;;
    -am-libname | --am-libname)
      prev=cfg_automake_libname
      ;;
    -am-noinst | --am-noinst)
      cfg_automake_noinst=true
      ;;
    -no-am-noinst | --no-am-noinst)
      cfg_automake_noinst=
      ;;

    -*)
      perror_exit "unrecognized option: ${opt}\nTry '$0 --help' for more information."
      ;;

    *=*)
      envvar=`expr "x${opt}" : 'x\([^=]*\)='`
      # Reject names that are not valid shell variable names.
      expr "x${envvar}" : ".*[^_${cr_alnum}]" >/dev/null && { echo "${me}: ERROR: invalid variable name: ${envvar}" >&2 ; exit 1 ; }
      optarg=`echo "${optarg}" | sed "s/'/'\\\\\\\\''/g"`
      eval "${envvar}='${optarg}'"
      export ${envvar}
      ;;

    *)
      perror_exit "this is not a valid option or argument: '${opt}'\nTry '$0 --help' for more information."
      ;;
  esac
  
  fi

  # if the previous option needs an argument, assign it
  if [ -n "${prev}" -a -n "${opt_isset}" ] ; then
    [ -n "${append}" ] && eval "append=\${${prev}}"
    if [ -n "${append}" ] ; then
      eval "${prev}=\"${append} \${opt}\""
    else
      eval "${prev}=\${opt}"
    fi
    prev=
    append=
  fi

done

[ -n "${prev}" ] && perror_exit "missing argument to --${prev}"


pprogress ""
pprogress "  SL - (parallel) sorting library (v${VERSION})"
pprogress ""

pdebug "options: $*"

# complete settings

# find-and-replace support
[ -z "${cfg_far}" ] && cfg_far="perl"

# config defaults (if empty try 'config/' and 'config.h')
if [ -z "${cfg_config}" ] ; then
  cfg_config="config/"
  [ ! -e "${cfg_config}" ] && cfg_config="config.h"
  [ ! -e "${cfg_config}" ] && cfg_config=
fi

# cfg_config
[ -d "${cfg_config}" ] && canonize_dir cfg_config

# cfg_config_not
if [ -n "${cfg_config_not}" ] ; then
  IFS_SAVE="${IFS}"
  IFS=", "
  set ${cfg_config_not}
  IFS="${IFS_SAVE}"
  cfg_config_not=" $* "
fi

# src_sl directories
canonize_dir src_sl src_sl_src src_sl_extra src_sl_scripts src_sl_adds

# dst_sl directories
canonize_dir dst_sl dst_sl_src dst_sl_extra dst_sl_scripts dst_sl_adds dst_if

[ -n "${cfg_source_ref_set}" ] && canonize_dir cfg_source_ref

# verify settings

# cfg_config
[ -n "${cfg_config}" -a ! -e "${cfg_config}" ] && perror_exit "configuration directory/file '${cfg_config}' does not exist!"

# src_sl
[ -n "${src_sl}" -a ! -d "${src_sl}" ] && perror_exit "source directory '${src_sl}' does not exist!"


# di_versions
if [ -n "${cfg_di_versions}" ] ; then
  case ${cfg_di_versions} in
    0|no|false)
      cfg_di_versions=
      ;;
    config|makefile)
      ;;
    *)
      perror_exit "invalid argument '${cfg_di_versions}' for option 'di-versions'"
  esac
fi

# sl_use_mpi
if [ -n "${cfg_makefile_sl_use_mpi}" ] ; then
  case ${cfg_makefile_sl_use_mpi} in
    0|no|false|ignore)
      cfg_makefile_sl_use_mpi=
      ;;
    1|yes|true|force)
      cfg_makefile_sl_use_mpi=true
      ;;
    *)
      perror_exit "invalid argument '${cfg_makefile_sl_use_mpi}' for option 'sl_use_mpi'"
  esac
fi

# cfg_local
[ -n "${cfg_local}" ] && cfg_source=ref

# cfg_source
[ -z "${cfg_source}" ] && cfg_source=ref
[ "${cfg_source}" = "sep" ] && cfg_source=separate

# ref settings
if [ -n "${cfg_source_ref_set}" ] ; then
  ref_sl=${cfg_source_ref}
  ref_sl_src=${ref_sl}
  [ -n "${src_sl_src_sub}" ] && ref_sl_src=${ref_sl_src}${src_sl_src}
  ref_sl_extra=${ref_sl}
  [ -n "${src_sl_extra_sub}" ] && ref_sl_extra=${ref_sl_extra}${src_sl_extra}
  ref_sl_scripts=${ref_sl}
  [ -n "${src_sl_scripts_sub}" ] && ref_sl_scripts=${ref_sl_scripts}${src_sl_scripts}
  ref_sl_adds=${ref_sl}
  [ -n "${src_sl_adds_sub}" ] && ref_sl_adds=${ref_sl_adds}${src_sl_adds}
fi


# verbose settings
print_settings pverbose " "


# build source directories
[ -n "${src_sl_src_sub}" ] && src_sl_src=${src_sl}${src_sl_src}
[ -n "${src_sl_extra_sub}" ] && src_sl_extra=${src_sl}${src_sl_extra}
[ -n "${src_sl_scripts_sub}" ] && src_sl_scripts=${src_sl}${src_sl_scripts}
[ -n "${src_sl_adds_sub}" ] && src_sl_adds=${src_sl}${src_sl_adds}


# complete commands
[ -z "${cmd_cat}" ] && cmd_cat=`which cat`
#[ -z "${cmd_cpp}" ] && cmd_cpp=`which cpp`
[ -z "${cmd_grep}" ] && cmd_grep=`which grep`
[ -z "${cmd_sed}" ] && cmd_sed=`which sed`
[ -z "${cmd_cp}" ] && cmd_cp=`which cp`
[ -z "${cmd_mkdir}" ] && cmd_mkdir=`which mkdir`

find_and_replace="${src_sl_scripts}find_and_replace.sh"
cmd_far="${find_and_replace} ${cfg_far}"


# verify commands
${cmd_sed} -e "s/0/0/" </dev/null >/dev/null 2>/dev/null || perror_exit "invalid sed command: '${cmd_sed}'"

${cmd_cat} </dev/null >/dev/null 2>/dev/null || perror_exit "invalid cat command: '${cmd_cat}'"

echo "X" | ${cmd_grep} "X" >/dev/null 2>/dev/null || perror_exit "invalid grep command: '${cmd_grep}'"

#${cmd_cpp} </dev/null >/dev/null 2>/dev/null || perror_exit "invalid cpp command: '${cmd_cpp}'"

[ -n "${cmd_cp}" ] || perror_exit "invalid cp command: '${cmd_cp}'"

[ -n "${cmd_mkdir}" ] || perror_exit "invalid mkdir command: '${cmd_mkdir}'"


# verbose commands
print_commands pverbose " "

# setup configuration
sub_config=config/
sub_include=include/
sub_base=base/
sub_base_mpi=base_mpi/
sub_tune=tune/
sub_tune_mpi=tune_mpi/

src_config=${src_sl_src}${sub_config}
src_include=${src_sl_src}${sub_include}
src_base=${src_sl_src}${sub_base}
src_base_mpi=${src_sl_src}${sub_base_mpi}
src_tune=${src_sl_src}${sub_tune}
src_tune_mpi=${src_sl_src}${sub_tune_mpi}

src_adds=${src_sl_adds}

src_config_env=${src_config}environment.h
src_config_tune=${src_config}tune.h

config_file=${cfg_config##*/}
config_dir=${cfg_config%${config_file}}

config_env=${config_dir}environment.h
config_tune=${config_dir}tune.h


# verbose: configuration
pverbose ""
pverbose " configuration:"
pverbose "  source environment file: ${src_config_env}"
pverbose "  source tune file: ${src_config_tune}"
pverbose "  config environment file: ${config_env}"
pverbose "  config tune file: ${config_tune}"

# config directory
if [ -d "${cfg_config}" ] ; then
  pverbose "  located configuration file(s):"

  # enumerate *.h files
  for file in ${cfg_config}*.h ; do
    # skip environment.h, tune.h and _tune.h files
    [ ! -f "${file}" -o "${file##*/}" = "environment.h" -o "${file##*/}" = "tune.h" -o "${file}" != "${file%_tune.h}" ] && continue

    all_config_files="${all_config_files}${file} "
    pverbose "   ${file}"
  done
else
  all_config_files=${cfg_config}
  pverbose "  config file: ${cfg_config}"
fi

#[ -z "${all_config_files}" ] && perror_exit "no configuration file found!"


# exit if 'testing'
[ -n "${cfg_testing}" ] && exit 0


# remove a local configuration
if [ -n "${cfg_clean}" ] ; then

  dst="${src_include}"

  pprogress ""
  pprogress -n " removing local files ..."

  rm -rf ${dst}sl_config.h ${dst}sl_environment.h ${dst}sl_tune.h ${dst}sl_tune_auto.h ${dst}sl_extra.h
  rm -rf autom4te.cache config.log config.status configure

  pprogress " done"
  pprogress ""
  pprogress ""

  exit 0
fi


create_sl_config()
{
  local dst_file="$1"
  local config_file="$2"
  
  pprogress -n "    ${dst_file} ..."

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"             > ${dst_file}
  else
    printf ""                             > ${dst_file}
  fi

  printf "\n"                             >> ${dst_file}
  printf "#ifndef __SL_CONFIG_H__\n"      >> ${dst_file}
  printf "#define __SL_CONFIG_H__\n"      >> ${dst_file}
  printf "\n"                             >> ${dst_file}
  ${cmd_cat} "${config_file}"             >> ${dst_file}
  printf "\n"                             >> ${dst_file}
  printf "#endif /* __SL_CONFIG_H__ */\n" >> ${dst_file}

  pprogress " done"
}

create_sl_environment()
{
  local dst_file="$1"
  local config_env="$2"
  local src_config_env="$3"

  pprogress -n "    ${dst_file} ..."

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"                  > ${dst_file}
  else
    printf ""                                  > ${dst_file}
  fi

  printf "\n"                                  >> ${dst_file}
  printf "#ifndef __SL_ENVIRONMENT_H__\n"      >> ${dst_file}
  printf "#define __SL_ENVIRONMENT_H__\n"      >> ${dst_file}
  printf "\n"                                  >> ${dst_file}
  # ... from $config_env/environment.h
  if [ -f "${config_env}" ] ; then
    ${cmd_cat} "${config_env}"                 >> ${dst_file}
    printf "\n"                                >> ${dst_file}
  fi
  # ... from $src_sl/config/environment.h
  if [ -f "${src_config_env}" ] ; then
    ${cmd_cat} "${src_config_env}"             >> ${dst_file}
    printf "\n"                                >> ${dst_file}
  fi
  printf "#endif /* __SL_ENVIRONMENT_H__ */\n" >> ${dst_file}

  pprogress " done"
}

create_sl_tune()
{
  local dst_file="$1"
  local cur_tune_file="$2"
  local config_tune_file="$3"
  local src_tune_file="$4"

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"           > ${dst_file}
  else
    printf ""                           > ${dst_file}
  fi
  pprogress -n "    ${dst_file} ..."

  printf "\n"                           >> ${dst_file}
  printf "#ifndef __SL_TUNE_H__\n"      >> ${dst_file}
  printf "#define __SL_TUNE_H__\n"      >> ${dst_file}
  printf "\n"                           >> ${dst_file}
  # ... from $config_dir/<prefix>_tune.h
  if [ -f "${cur_tune_file}" ] ; then
    cat "${cur_tune_file}"              >> ${dst_file}
    printf "\n"                         >> ${dst_file}
  fi
  # ... from $config_dir/tune.h
  if [ -f "${config_tune_file}" ] ; then
    cat ${config_tune_file}             >> ${dst_file}
    printf "\n"                         >> ${dst_file}
  fi
  # ... from $src_sl/config/tune.h
  if [ -f "${src_tune_file}" ] ; then
    cat ${src_tune_file}                >> ${dst_file}
    printf "\n"                         >> ${dst_file}
  fi
  printf "#endif /* __SL_TUNE_H__ */\n" >> ${dst_file}

  pprogress " done"
}

create_sl_tune_auto()
{
  local dst_file="$1"
  
  pprogress -n "    ${dst_file} ..."

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"                > ${dst_file}
  else
    printf ""                                > ${dst_file}
  fi

  printf "\n"                                >> ${dst_file}
  printf "#ifndef __SL_TUNE_AUTO_H__\n"      >> ${dst_file}
  printf "#define __SL_TUNE_AUTO_H__\n"      >> ${dst_file}
  printf "\n"                                >> ${dst_file}
  printf "\n"                                >> ${dst_file}
  printf "#endif /* __SL_TUNE_AUTO_H__ */\n" >> ${dst_file}

  pprogress " done"
}

create_sl_extra()
{
  local dst_file="$1"
  local x
  
  pprogress -n "    ${dst_file} ..."

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"            > ${dst_file}
  else
    printf ""                            > ${dst_file}
  fi

  printf "\n"                            >> ${dst_file}
  printf "#ifndef __SL_EXTRA_H__\n"      >> ${dst_file}
  printf "#define __SL_EXTRA_H__\n"      >> ${dst_file}
  printf "\n"                            >> ${dst_file}
  if [ -n "${cfg_extra_have_h}" ] ; then
    printf "\n"                          >> ${dst_file}
    printf "#ifdef SL_EXTRA_PREFIX\n"    >> ${dst_file}
    printf "# define ZMPI_PREFIX  SL_EXTRA_PREFIX\n" >> ${dst_file}
    printf "#endif\n"                    >> ${dst_file}
    printf "\n"                          >> ${dst_file}
    printf "\n"                          >> ${dst_file}
    for x in ${cfg_extra_have_h} ; do
      printf "#define HAVE_`echo ${x} | tr [a-z] [A-Z]`_H\n" >> ${dst_file}
    done
    printf "\n"                          >> ${dst_file}
  fi
  printf "\n"                            >> ${dst_file}
  printf "#endif /* __SL_EXTRA_H__ */\n" >> ${dst_file}

  pprogress " done"
}


create_file_far_sed()
{
  local src_file="$1"
  local dst_file="$2"
  local far_script="$3"
  local sed_cmd="$4"
  
  [ -z "${dst_file##*/}" ] && dst_file="${dst_file}${src_file##*/}"

  pprogress -n "    ${dst_file} ..."

  if [ -n "${far_script}" ] ; then
    if [ -n "${sed_cmd}" ] ; then
      ${cmd_far} "exec" "${far_script}" "${src_file}" | ${cmd_sed} -e "${sed_cmd}" > ${dst_file}
    else
      ${cmd_far} "exec" "${far_script}" "${src_file}" > ${dst_file}
    fi
  else
    if [ -n "${sed_cmd}" ] ; then
      ${cmd_sed} -e "${sed_cmd}" "${src_file}" > ${dst_file}
    else
      ${cmd_cat} "${src_file}" > "${dst_file}"
    fi
  fi

  pprogress " done"
}

cat_files()
{
  local files="$*"

  for f in ${files} ; do
    [ -f "${f}" ] && ${cmd_cat} "${f}"
  done
}

include_include_files()
{
  local incdir="$1"
  shift
  local files="$*"

  for f in ${files} ; do
    sed_script="${sed_script} -e '/#include \"${f}\"/r ${incdir}${f}' -e '/#include \"${f}\"/d'"
  done

  eval ${cmd_sed} ${sed_script}
}

create_sl_interface()
{
  local dst_file=$1
  local cur_config_name=$2
  local cur_config_file=$3
  local cur_tune_file=$4
  local config_tune_file=$5
  local src_tune_file=$6
  local src_incdir=$7
  local far_script=$8

  local uppercase_config=`printf "${cur_config_name}" | tr a-z A-Z`

  if [ -f "${cfg_license_file}" ] ; then
    cat "${cfg_license_file}"                     > ${dst_file}
  else
    printf ""                                     > ${dst_file}
  fi

  printf "\n"                                     >> ${dst_file}
  printf "#ifndef __SL_${uppercase_config}_H__\n" >> ${dst_file}
  printf "#define __SL_${uppercase_config}_H__\n" >> ${dst_file}
  printf "\n"                                     >> ${dst_file}
  printf "#ifdef SL_USE_MPI\n"                    >> ${dst_file}
  printf " #include <mpi.h>\n"                    >> ${dst_file}
  printf "#endif /* SL_USE_MPI */\n"              >> ${dst_file}
  printf "\n"                                     >> ${dst_file}
  printf "#define SL_PROTO(_f_)  _f_"             >> ${dst_file}
  printf "\n"                                     >> ${dst_file}
  cat_files "${cur_config_file}" \
            "${cur_tune_file}" "${config_tune_file}" "${src_tune_file}" \
            "${src_incdir}/sl_config_global.h" \
            "${src_incdir}/sl_config_intern.h" \
            "${src_incdir}/sl_tune_intern.h" \
            "${src_incdir}/sl_rti_tids.h" \
            "${src_incdir}/spec_public_conf.h" \
            "${src_incdir}/spec_public.h" \
            "${src_incdir}/sl_types.h" \
            "${src_incdir}/sl_adds.h" \
            "${src_incdir}/sl_globals.h" \
            "${src_incdir}/sl_protos.h" \
   | ${cmd_sed} -e '$a\
#ifdef SL_USE_MPI\
' \
   | ${cmd_cat} - "${src_incdir}/sl_protos_mpi.h" \
   | ${cmd_sed} -e '$a\
#endif \/\* SL_USE_MPI \*\/\
' \
   | ${cmd_grep} -v "__" \
   | include_include_files "${src_incdir}/" "sl_context_struct.h" \
   | ${cmd_far} exec "${far_script}" \
   | ${cmd_sed} -e '/\/\*$/,/ \*\/$/d'                 >> ${dst_file}
  printf "\n"                                          >> ${dst_file}
  printf "#undef SL_PROTO"                             >> ${dst_file}
  printf "\n"                                          >> ${dst_file}
  printf "#endif /* __SL_${uppercase_config}_H__ */\n" >> ${dst_file}
}

create_includes()
{
  local src=$1
  local dst=$2
  local sub=$3
  local file src_file dst_file

  mkdir_save "${dst}${sub}"

  for src_file in `echo ${src}${sub}*.h | tr ' ' '\n' | sort` ; do
    [ ! -f ${src_file} ] && continue

    file="${src_file##*/}"
    dst_file="${dst}${sub}${file}"

    [ "${file}" = "sl_config.h" ] && continue
    [ "${file}" = "sl_environment.h" ] && continue
    [ "${file}" = "sl_tune.h" ] && continue
    [ "${file}" = "sl_auto_tune.h" ] && continue
    [ "${file}" = "sl_extra.h" ] && continue

    pdebug "file: ${file}"
    pdebug "src_file: ${src_file}"
    pdebug "dst_file: ${dst_file}"

    create_file_far_sed "${src_file}" "${dst_file}"
  done
}

create_sources()
{
  local src=$1
  local dst=$2
  local sub=$3
  local pat=$4
  local prefix=$5
  local file src_file dst_file

  mkdir_save "${dst}${sub}"

  for src_file in `echo ${src}${sub}${pat} | tr ' ' '\n' | sort` ; do
    [ ! -f ${src_file} ] && continue

    file="${src_file##*/}"
    dst_file="${dst}${sub}${prefix}${file}"

    pdebug "file: ${file}"
    pdebug "src_file: ${src_file}"
    pdebug "dst_file: ${dst_file}"

    create_file_far_sed "${src_file}" "${dst_file}"

#    if [ "${cfg_di_versions}" = "config" ] ; then
#      # base file di-version, rename variable & function names!
#      create_file_far_sed "${file}" "${dst_file_di}" "${far_script_varfuncs_di}" "1i#define SL_DATA_IGNORE"
#    fi
  done
}

create_scripts()
{
  local src=$1
  local dst=$2

  local dst_scripts="${dst}${dst_sl_scripts}"

  local file src_file dst_file

#  pdebug "src_sl_scripts: ${src_sl_scripts}"
#  pdebug "dst_scripts: ${dst_scripts}"

  mkdir_save "${dst_scripts}"

  for src_file in ${src_sl_scripts}* ; do
    [ ! -f ${src_file} ] && continue

    file="${src_file##*/}"
    dst_file="${dst_scripts}${file}"

    pdebug "file: ${file}"
    pdebug "src_file: ${src_file}"
    pdebug "dst_file: ${dst_file}"

    create_file_far_sed "${src_file}" "${dst_file}"

    chmod u+x "${dst_file}"
  done
}

create_adds()
{
  local src=$1
  local dst=$2

  local dst_adds="${dst}${dst_sl_adds}"

  pdebug "src_sl_adds: ${src_sl_adds}"
  pdebug "dst_adds: ${dst_adds}"

  mkdir_save "${dst_adds}"

  create_file_far_sed "${src_sl_adds}names_oxl_all_functions" "${dst_adds}"
}


# source: ref
#if [ "${cfg_source}" = "ref" ] ; then
#fi

mkdir_save "${dst_sl}"

# source: single
if [ "${cfg_source}" = "single" ] ; then

  pprogress ""
  pprogress " installing single source ..."

  current_dst="${dst_sl}"

  mkdir_save "${current_dst}"

  if [ -n "${dst_sl_src_sub}" ] ; then
    current_dst_sl_src="${current_dst}${dst_sl_src}"
  else
    current_dst_sl_src="${dst_sl_src}"
  fi
  
  pdebug "current_dst_sl_src: ${current_dst_sl_src}"

  # include files
  pprogress ""
  pprogress "   creating source headers ..."

  create_includes "${src_sl_src}" "${current_dst_sl_src}" "${sub_include}"

  # base files
  pprogress ""
  pprogress "   creating base source files ..."

  create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_base}" "*.c"
  create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_base_mpi}" "*.c"

  # tuning
  if [ -n "${cfg_tuning}" ] ; then
    # tune files
    pprogress ""
    pprogress "   creating tune source files ..."

    create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_tune}" "*.[ch]"
    create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_tune_mpi}" "*.[ch]"

    pprogress ""
    pprogress "   creating automatical tuning support ..."

    # sl_tune.sh
    create_file_far_sed "${src_sl}sl_tune.sh" "${current_dst}"
  fi

  if [ -n "${dst_sl_extra_sub}" ] ; then
    current_dst_sl_extra="${current_dst}${dst_sl_extra}"
  else
    current_dst_sl_extra="${dst_sl_extra}"
  fi
  
  pdebug "current_dst_sl_extra: ${current_dst_sl_extra}"

  # extra
  if [ -n "${cfg_extra}" ] ; then
    pprogress ""
    pprogress "   creating extra files ..."

    create_sources "${src_sl_extra}" "${current_dst_sl_extra}" "" "*"

    for x in ${cfg_extras} include ; do
      create_sources "${src_sl_extra}" "${current_dst_sl_extra}" "${x}/" "*.[hc]"
    done
  fi

  # scripts
  if [ -n "${cfg_scripts}" ] ; then
    pprogress ""
    pprogress "   creating script files ..."

    create_scripts "DUMMY" "${current_dst}"
    create_adds "DUMMY" "${current_dst}"
  fi
fi


config_names=
config_not_names=

pprogress ""
pprogress " installing configurations ..."

for current_config_file in ${all_config_files} ; do

  current_config_filename="${current_config_file##*/}"
  current_config_name="${current_config_filename%.h}"

  case ${cfg_config_not} in
    *\ ${current_config_name}\ *)
      config_names_not="${config_names_not} ${current_config_name}"
      ;;
    *)
      config_names="${config_names} ${current_config_name}"
      ;;
  esac

  current_tune_file="${config_dir}${current_config_name}_tune.h"

  if [ -n "${cfg_local}" ] ; then
    current_dst="${dst_sl}"
  else
    current_dst="${dst_sl}${dstdir_prefix}${current_config_name}/"
  fi

  if [ -n "${cfg_source_rename}" ] ; then
    current_prefix="${current_config_name}_"
  else
    current_prefix=
  fi

  pprogress ""
  pprogress "  installing '${current_config_name}' ..."
  pverbose  "   config file: ${current_config_file}"
  pverbose  "   tune file: ${current_tune_file}"
  pverbose  "   destination directory: ${current_dst}"

  mkdir_save "${current_dst}"

  if [ -n "${dst_sl_src_sub}" ] ; then
    current_dst_sl_src="${current_dst}${dst_sl_src}"
  else
    if [ "${cfg_source}" = "separate" ] ; then
      perror "dst_sl_src_sub cannot be 'not true'!"
    else
      current_dst_sl_src="${dst_sl_src}"
    fi
  fi

  pdebug "current_dst_sl_src: ${current_dst_sl_src}"

  if [ -n "${cfg_local}" ] ; then
    current_dst_cfg="${current_dst_sl_src}${sub_include}"
  else
    if [ "${cfg_source}" = "separate" ] ; then
      current_dst_cfg="${current_dst_sl_src}${sub_include}"
    else
      current_dst_cfg="${current_dst}"
    fi
  fi

  pdebug "current_dst_cfg: ${current_dst_cfg}"

  pprogress ""
  pprogress "   creating source headers ..."

  mkdir_save "${current_dst_cfg}"

  # sl_config.h
  create_sl_config "${current_dst_cfg}sl_config.h" "${current_config_file}"

  # sl_environment.h
  create_sl_environment "${current_dst_cfg}sl_environment.h" "${config_env}" "${src_config_env}"

  # sl_tune.h
  create_sl_tune "${current_dst_cfg}sl_tune.h" "${current_tune_file}" "${config_tune}" "${src_config_tune}"

  # sl_tune_auto.h
  if [ -n "${cfg_tune_auto_dummy}" ] ; then
    create_sl_tune_auto "${current_dst_cfg}sl_tune_auto.h"
  fi

  # sl_extra.h
  if [ -n "${cfg_extra}" ] ; then
    create_sl_extra "${current_dst_cfg}sl_extra.h"
  fi

  if [ "${cfg_source}" = "separate" ] ; then

    create_includes "${src_sl_src}" "${current_dst_sl_src}" "${sub_include}"

    # base files
    pprogress ""
    pprogress "   creating base source files ..."

    create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_base}" "*.c" "${current_prefix}"
    create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_base_mpi}" "*.c" "${current_prefix}"

    # tuning
    if [ -n "${cfg_tuning}" ] ; then
      # tune files
      pprogress ""
      pprogress "   creating tune source files ..."

      create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_tune}" "*.[ch]" "${current_prefix}"
      create_sources "${src_sl_src}" "${current_dst_sl_src}" "${sub_tune_mpi}" "*.[ch]" "${current_prefix}"

      pprogress ""
      pprogress "   creating automatical tuning support ..."

      # sl_tune.sh
      create_file_far_sed "${src_sl}sl_tune.sh" "${current_dst}"
    fi

    # scripts
    if [ -n "${cfg_scripts}" ] ; then
      pprogress ""
      pprogress "   creating script files ..."

      create_scripts "DUMMY" "${current_dst}"
      create_adds "DUMMY" "${current_dst}"
    fi
  fi

  # interface
  if [ -n "${cfg_interface}" ] ; then

    pprogress ""
    pprogress "   creating library interface ..."

    dst_include="${dst_sl}"

    if [ -n "${dst_if_sub}" ] ; then
      dst_include="${dst_include}${dst_if}"
    else
      dst_include="${dst_if}"
    fi

    pdebug "dst_include: ${dst_include}"

    mkdir_save "${dst_include}"

    # sl_<prefix>.h
    dst_interface="${dst_include}${header_prefix}${current_config_name}.h"
    pprogress -n "    ${dst_interface} ..."

    far_script_interface="${current_dst}far_script_interface"

    ${cmd_far} begin                                                                           > ${far_script_interface}

    # make rename-script for type to <prefix>_type
    ${cmd_far} mod_list 1 "${src_adds}names_types" ":" "${cfg_prefix}${current_config_name}_:"             >> ${far_script_interface}
    ${cmd_far} mod_list 1 "${src_adds}spec_names_types" ":" "${cfg_prefix}${current_config_name}_:"        >> ${far_script_interface}
    # make rename-script for macro to <prefix>_macro
    ${cmd_far} mod_list 1 "${src_adds}names_macros" ":" "${cfg_prefix}${current_config_name}_:"            >> ${far_script_interface}
    ${cmd_far} mod_list 1 "${src_adds}spec_names_macros" ":" "${cfg_prefix}${current_config_name}_:"       >> ${far_script_interface}
    # make rename-script for [variable|function] to <prefix>_[variable|function]
    ${cmd_far} mod_list 1 "${src_adds}names_oxl_all_variables" ":" "${cfg_prefix}${current_config_name}_:" >> ${far_script_interface}
    ${cmd_far} mod_list 1 "${src_adds}names_oxl_all_functions" ":" "${cfg_prefix}${current_config_name}_:" >> ${far_script_interface}
    # make rename-script for QWERTZ_macro to macro
    ${cmd_far} mod_list 1 "${src_adds}names_macros" "QWERTZ_:" ":"                                         >> ${far_script_interface}

    ${cmd_far} end                                                                                         >> ${far_script_interface}

    create_sl_interface "${dst_interface}" "${current_config_name}" "${current_config_file}" "${current_tune_file}" "${config_tune}" "${src_config_tune}" "${src_include}" "${far_script_interface}"

    [ -z "${cfg_debug}" ] && rm -f "${far_script_interface}"

    pprogress " done"

  fi

  pprogress ""
  pprogress "  installing '${current_config_name}' in '${dst_current}' finished!"
  
  [ -n "${cfg_local}" ] && break
  
done

pprogress ""
pprogress " installing configurations in '${dst_sl}' finished!"


# setup makefile-find-and-replace string
makefile_makefile_in="${cfg_makefile_makefile_in}"


# Makefile
if [ -n "${cfg_makefile}" ] ; then
  pprogress ""
  pprogress " creating Makefile ..." >&3

  makefile_sl_use_mpi="${cfg_makefile_sl_use_mpi}"
  makefile_di_versions="${cfg_di_versions}"
  makefile_far="${cfg_far}"
  makefile_target="${cfg_makefile_target}"
  makefile_quiet="${cfg_makefile_quiet}"
  makefile_bulk_ar="${cfg_makefile_bulk_ar}"
  makefile_prefix="${cfg_prefix}"
  makefile_extra="${cfg_extra}"
  makefile_extra_prefix="${cfg_extra_prefix}"
  if [ "${cfg_source}" = "ref" ] ; then
    if [ -n "${cfg_source_ref_set}" ] ; then
      makefile_srcbase="${ref_sl}"
    else
      makefile_srcbase="${src_sl}"
    fi
  fi
  makefile_all="${cfg_makefile_all}"

  makefile_wrapper_src="${cfg_makefile_wrapper_src}"
  makefile_wrapper_prefix="${cfg_makefile_wrapper_prefix}"
  makefile_wrapper_mpi_src="${cfg_makefile_wrapper_mpi_src}"
  makefile_wrapper_mpi_prefix="${cfg_makefile_wrapper_mpi_prefix}"

  makefile_exec="${cfg_makefile_exec}"

  makefile_exec_src="${cfg_makefile_exec_src}"
  makefile_exec_prefix="${cfg_makefile_exec_prefix}"
  makefile_exec_mpi_src="${cfg_makefile_exec_mpi_src}"
  makefile_exec_mpi_prefix="${cfg_makefile_exec_mpi_prefix}"

  makefile_incdir="${cfg_makefile_incdir}"
  makefile_libdir="${cfg_makefile_libdir}"

  create_far_entry()
  {
    pverbose "    setting '$1=$2'"
    ${cmd_far} single 0 "#$1#" "$2" >> $3
  }

  far_script_makefile="${dst}far_makefile"
  ${cmd_far} begin > ${far_script_makefile}

  create_far_entry "MAKEFILE_IN" "${makefile_makefile_in}" "${far_script_makefile}"
  create_far_entry "SL_USE_MPI" "${makefile_sl_use_mpi}" "${far_script_makefile}"
#  create_far_entry "MAX_NDATA" "${makefile_max_ndata}" "${far_script_makefile}"
  create_far_entry "DI_VERSIONS" "${makefile_di_versions}" "${far_script_makefile}"
  create_far_entry "FAR" "${makefile_far}" "${far_script_makefile}"
  create_far_entry "TARGET" "${makefile_target}" "${far_script_makefile}"
  create_far_entry "QUIET" "${makefile_quiet}" "${far_script_makefile}"
  create_far_entry "BULK_AR" "${makefile_bulk_ar}" "${far_script_makefile}"
  create_far_entry "SRCBASE" "${makefile_srcbase}" "${far_script_makefile}"
  create_far_entry "ALL" "${makefile_all}" "${far_script_makefile}"

  create_far_entry "PREFIX" "${makefile_extra}" "${far_script_makefile}"

  create_far_entry "EXTRA" "${makefile_extra}" "${far_script_makefile}"
  create_far_entry "EXTRA_PREFIX" "${makefile_extra_prefix}" "${far_script_makefile}"

  create_far_entry "WRAPPER_SRC" "${makefile_wrapper_src}" "${far_script_makefile}"
  create_far_entry "WRAPPER_PREFIX" "${makefile_wrapper_prefix}" "${far_script_makefile}"
  create_far_entry "WRAPPER_MPI_SRC" "${makefile_wrapper_mpi_src}" "${far_script_makefile}"
  create_far_entry "WRAPPER_MPI_PREFIX" "${makefile_wrapper_mpi_prefix}" "${far_script_makefile}"

  create_far_entry "EXEC" "${makefile_exec}" "${far_script_makefile}"

  create_far_entry "EXEC_SRC" "${makefile_exec_src}" "${far_script_makefile}"
  create_far_entry "EXEC_PREFIX" "${makefile_exec_prefix}" "${far_script_makefile}"
  create_far_entry "EXEC_MPI_SRC" "${makefile_exec_mpi_src}" "${far_script_makefile}"
  create_far_entry "EXEC_MPI_PREFIX" "${makefile_exec_mpi_prefix}" "${far_script_makefile}"

  create_far_entry "INCDIR" "${makefile_incdir}" "${far_script_makefile}"
  create_far_entry "LIBDIR" "${makefile_libdir}" "${far_script_makefile}"

  ${cmd_far} end >> ${far_script_makefile} 

  create_file_far_sed "${src_sl}Makefile" "${dst_sl}Makefile" "${far_script_makefile}"

  if [ -n "${cfg_makefile_fixed}" ] ; then
    ${cmd_sed} -i -e "s/\$(foreach.*MODULES.*$/FIXED_REPLACE/" "${dst_sl}Makefile"
    ${cmd_sed} -i -e "s/^FIXED_REPLACE/ifneq (\$(MODULES),)\nFIXED_REPLACE/" "${dst_sl}Makefile"
    for c in ${config_names} ; do
      ${cmd_sed} -i -e "s/^FIXED_REPLACE/ ifeq (\$(findstring sl_${c},\$(\$(mod)_MODULES_NOT)),)\nFIXED_REPLACE/" "${dst_sl}Makefile"
      ${cmd_sed} -i -e "s/^FIXED_REPLACE/  mod:=sl_${c}\/\nFIXED_REPLACE/" "${dst_sl}Makefile"
      ${cmd_sed} -i -e "s/^FIXED_REPLACE/  include Makefile\nFIXED_REPLACE/" "${dst_sl}Makefile"
      ${cmd_sed} -i -e "s/^FIXED_REPLACE/  mod:=\$(firstword \$(mod_list))\nFIXED_REPLACE/" "${dst_sl}Makefile"
      ${cmd_sed} -i -e "s/^FIXED_REPLACE/ endif\nFIXED_REPLACE/" "${dst_sl}Makefile"
    done
    ${cmd_sed} -i -e "s/^FIXED_REPLACE/endif\n/" "${dst_sl}Makefile"
  fi

  [ -z "${cfg_debug}" ] && rm -f ${far_script_makefile}

fi


# Makefile.in
if [ -n "${cfg_makefile_in}" ] ; then
  mf_in_src="${src_sl}${cfg_makefile_in}"
  mf_in_dst="${dst_sl}${makefile_makefile_in}"
  if [ -f "${mf_in_src}" ] ; then
    pprogress ""
    pprogress " creating ${mf_in_dst} ..."
    create_file_far_sed "${mf_in_src}" "${mf_in_dst}"
  fi
fi


# Autoconf files
if [ -n "${cfg_autoconf}" ] ; then
  pprogress ""
  pprogress " creating Autoconf files ..."
  create_file_far_sed "${src_sl}Makefile.in.in" "${dst_sl}"
  create_file_far_sed "${src_sl}configure.ac" "${dst_sl}"
  ${cmd_mkdir} "${dst_sl}build-aux"
  create_file_far_sed "${src_sl}build-aux/acx_mpi.m4" "${dst_sl}build-aux/"
  create_file_far_sed "${src_sl}build-aux/install-sh" "${dst_sl}build-aux/"
#  [ -f "${src_sl}configure" ] && create_file_far_sed "${src_sl}configure" "${dst_sl}"
fi


# Automake file
if [ -n "${cfg_automake}" ] ; then

  pprogress ""
  pprogress " creating Automake file(s) ..."

  pprogress "    ${dst_sl}Makefile.am ..."

  rm -f "${dst_sl}Makefile.am"

  exec >> "${dst_sl}Makefile.am"

  echo ""
  echo "lib_LIBRARIES ="
  echo "noinst_LIBRARIES ="
  echo ""
  if [ -n "${cfg_automake_noinst}" ] ; then
    echo -n "noinst_"
  else
    echo -n "lib_"
  fi
  echo "LIBRARIES += ${cfg_automake_libname}"
  echo ""

  if [ "${cfg_source}" = "ref" ] ; then
    current_dst_sl_src="${dst_sl}"
    sub_src=${src_sl_src}
    current_dst_sl_extra="${dst_sl}"
    sub_extra=${src_sl_extra}
#    if [ -n "${cfg_source_ref_set}" ] ; then
#      sub_src=${cfg_source_ref}
#    else
#      sub_src=${src_sl_src}
#    fi
  elif [ "${cfg_source}" = "single" ] ; then
    if [ -n "${dst_sl_src_sub}" ] ; then
      current_dst_sl_src="${dst_sl}"
      sub_src="${dst_sl_src}"
    else
      current_dst_sl_src="${dst_sl_src}"
      sub_src=
    fi
    if [ -n "${dst_sl_extra_sub}" ] ; then
      current_dst_sl_extra="${dst_sl}"
      sub_extra="${dst_sl_extra}"
    else
      current_dst_sl_extra="${dst_sl_extra}"
      sub_extra=
    fi
  fi

  pdebug "current_dst_sl_src: ${current_dst_sl_src}"
  pdebug "current_dst_sl_extra: ${current_dst_sl_extra}"

  have_wrapper_sources_c=
  for f in ${dst_sl}*.c ; do
    [ ! -e "${f}" ] && continue
    have_wrapper_sources_c="yes"
  done

  echo -n "sl_sub_libs ="
  for c in ${config_names} ; do
    echo -n " libsl_${c}.a"
  done
  if [ -n "${cfg_extra}" ] ; then
    echo -n " libsl_extra.a"
  fi
  if [ -n "${have_wrapper_sources_c}" ] ; then
    echo -n " libsl_mpiwrap.a"
  fi
  echo ""
  echo ""
  echo "noinst_LIBRARIES += \$(sl_sub_libs)"
  echo ""
  echo "srcdir_sl=\$(srcdir)"
  if [ "${cfg_source}" = "ref" -o "${cfg_source}" = "single" ] ; then
    if [ -n "${cfg_source_ref_set}" ] ; then
      echo "srcdir_sl_source=\$(srcdir)/${ref_sl_src%%/}"
      echo "srcdir_sl_extra=\$(srcdir)/${ref_sl_extra%%/}"
    else
      echo "srcdir_sl_source=\$(srcdir)/${sub_src%%/}"
      echo "srcdir_sl_extra=\$(srcdir)/${sub_extra%%/}"
    fi
  fi
  echo ""
  for c in ${config_names} ; do
    echo "libsl_${c}_a_CPPFLAGS = \\"
    echo -n "  -DSL_PREFIX=${cfg_prefix}${c}_ -DSL_USE_MPI"
    if [ "${cfg_source}" = "ref" -o "${cfg_source}" = "single" ] ; then
      echo -n " -I\$(srcdir_sl_source)/include -I\$(srcdir_sl)/sl_${c}"
    elif [ "${cfg_source}" = "separate" ] ; then
      echo -n " -I\$(srcdir_sl)/sl_${c}/${dst_sl_src}${sub_include}"
    fi
    echo " -DSL_EXTRA -DSL_EXTRA_PREFIX=${cfg_extra_prefix} -I\$(srcdir_sl_extra)/include"
    echo ""
  done
  if [ -n "${cfg_extra}" ] ; then
    echo "libsl_extra_a_CPPFLAGS = \\"
    echo -n "  -DZ_PREFIX=${cfg_extra_prefix} -DZMPI_PREFIX=${cfg_extra_prefix}"
    for x in ${cfg_extra_include} ; do
      echo -n " -DHAVE_`echo ${x} | tr [a-z] [A-Z]`_H -I\$(srcdir_sl_extra)/${x}"
    done
    echo ""
    echo ""
  fi
  if [ -n "${have_wrapper_sources_c}" ] ; then
    echo "libsl_mpiwrap_a_CPPFLAGS = \\"
    echo -n "  -DSL_USE_MPI -I\$(srcdir_sl)/include"
    for c in ${config_names_not} ; do
      echo -n " -DNOT_sl_${c}"
    done
    if [ -n "${cfg_extra}" ] ; then
      echo -n " -DZMPI_PREFIX=${cfg_extra_prefix}"
      for x in ${cfg_extra_wrap_have_h} ; do
        echo -n " -DHAVE_`echo ${x} | tr [a-z] [A-Z]`_H"
      done
      echo -n " -I\$(srcdir_sl_extra)/include"
    fi
    echo ""
  fi
  echo ""
  echo `echo "${cfg_automake_libname}_SOURCES =" | tr -- -. _`
  echo ""
  echo "${cfg_automake_libname}: \$(sl_sub_libs)"
  echo "	rm -f \$@"
  echo "	\$(AR) \$(ARFLAGS) \$@ *.o"
  echo "	\$(RANLIB) \$@"
  if [ "${cfg_source}" = "ref" -o "${cfg_source}" = "single" ] ; then
    pdebug "current_dst_sl_src: ${current_dst_sl_src}"
    pdebug "sub_src: ${sub_src}"
    echo ""
    echo -n "sl_SOURCE ="
    for f in `echo ${current_dst_sl_src}${sub_src}${sub_base}*.c ${current_dst_sl_src}${sub_src}${sub_base_mpi}*.c ${current_dst_sl_src}${sub_src}${sub_include}*.h | tr ' ' '\n' | sort` ; do
      [ ! -e "${f}" ] && continue
      echo " \\"
      echo -n "  \$(srcdir_sl_source)/${f#${current_dst_sl_src}${sub_src}}"
    done
    if [ -n "${cfg_extra}" ] ; then
      for f in `echo ${current_dst_sl_extra}${sub_extra}include/*.h | tr ' ' '\n' | sort` ; do
        pdebug "file: ${f}"
        [ ! -e "${f}" ] && continue
        echo " \\"
        echo -n "  \$(srcdir_sl_extra)/${f#${current_dst_sl_extra}${sub_extra}}"
      done
    fi
    echo ""
  fi
  for c in ${config_names} ; do
    echo ""
    echo -n "libsl_${c}_a_SOURCES ="
    for f in `echo ${dst_sl}*.h | tr ' ' '\n' | sort` ; do
      [ ! -e "${f}" ] && continue
      echo " \\"
      echo -n "  \$(srcdir_sl)/${f#${dst_sl}}"
    done
    for f in `echo ${dst_sl}sl_${c}/*.h | tr ' ' '\n' | sort` ; do
      [ ! -e "${f}" ] && continue
      echo " \\"
      echo -n "  \$(srcdir_sl)/${f#${dst_sl}}"
    done
    if [ "${cfg_source}" = "ref" -o "${cfg_source}" = "single" ] ; then
      echo " \\"
      echo "  \$(sl_SOURCE)"
    else
      for f in `echo ${dst_sl}sl_${c}/${dst_sl_src}${sub_include}*.h ${dst_sl}sl_${c}/${dst_sl_src}${sub_base}*.c ${dst_sl}sl_${c}/${dst_sl_src}${sub_base_mpi}*.c | tr ' ' '\n' | sort` ; do
        [ ! -e "${f}" ] && continue
        echo " \\"
        echo -n "  \$(srcdir_sl)/${f#${dst_sl}}"
      done
    fi
  done
  if [ -n "${cfg_extra}" ] ; then
    pdebug "current_dst_sl_extra: ${current_dst_sl_extra}"
    pdebug "sub_extra: ${sub_extra}"
    echo ""
    echo -n "libsl_extra_a_SOURCES ="
    for x in ${cfg_extras} ; do
      for f in `echo ${current_dst_sl_extra}${sub_extra}${x}/*.[hc] | tr ' ' '\n' | sort` ; do
        pdebug "file: ${f}"
        [ ! -e "${f}" ] && continue
        echo " \\"
        echo -n "  \$(srcdir_sl_extra)/${f#${current_dst_sl_extra}${sub_extra}}"
      done
    done
    echo ""
  fi
  echo ""
  if [ -n "${have_wrapper_sources_c}" ] ; then
    echo -n "libsl_mpiwrap_a_SOURCES ="
    for f in `echo ${dst_sl}include/*.h ${dst_sl}*.[hc] | tr ' ' '\n' | sort`; do
      [ ! -e "${f}" ] && continue
      echo " \\"
      echo -n "  \$(srcdir_sl)/${f#${dst_sl}}"
    done
    echo ""
  else
    echo -n "noinst_HEADERS ="
    for f in `echo ${dst_sl}include/*.h ${dst_sl}*.h | tr ' ' '\n' | sort` ; do
      [ ! -e "${f}" ] && continue
      echo " \\"
      echo -n "  \$(srcdir_sl)/${f#${dst_sl}}"
    done
    echo ""
  fi

  exec >&-
fi

pprogress ""
pprogress ""

exit 0
