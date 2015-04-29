
lib_LTLIBRARIES =
noinst_LTLIBRARIES =

noinst_LTLIBRARIES += libsl_pepc.la

sl_sub_libs = libfcs_pepc_pepckeys.la libfcs_pepc_extra.la libfcs_pepc_mpiwrap.la

noinst_LTLIBRARIES += $(sl_sub_libs)

libfcs_pepc_pepckeys_la_CPPFLAGS = \
  -DSL_PREFIX=fcs_pepckeys_ -DSL_USE_MPI -I$(srcdir)/src/include -I$(srcdir)/sl_pepckeys -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_pepc_ -I$(srcdir)/extra/include

libfcs_pepc_extra_la_CPPFLAGS = \
  -DZ_PREFIX=fcs_pepc_ -DZMPI_PREFIX=fcs_pepc_ -DHAVE_Z_PACK_H -I$(srcdir)/extra/z_pack -DHAVE_ZMPI_LOCAL_H -I$(srcdir)/extra/zmpi_local -DHAVE_ZMPI_TOOLS_H -I$(srcdir)/extra/zmpi_tools

libfcs_pepc_mpiwrap_la_CPPFLAGS = \
  -DSL_USE_MPI -I$(srcdir)/include -DNOT_sl_pepcparts -DZMPI_PREFIX=fcs_pepc_ -DHAVE_ZMPI_LOCAL_H -DHAVE_ZMPI_TOOLS_H -DHAVE_ZMPI_ATAIP_H -DHAVE_ZMPI_ATASP_H -I$(srcdir)/extra/include

libsl_pepc_la_SOURCES =

libsl_pepc_la_LIBADD = $(sl_sub_libs)

sl_SOURCE = \
  src/base/base.c \
  src/base_mpi/base_mpi.c \
  src/base_mpi/dash.c \
  src/base_mpi/spec.c \
  src/base/z_crc32.c \
  src/base/z_pack.c \
  src/include/dash_aux_heap.h \
  src/include/dash_aux_static.h \
  src/include/dash_common.h \
  src/include/dash_conf.h \
  src/include/dash_core.h \
  src/include/dash_exec.h \
  src/include/dash_exec_sl.h \
  src/include/dash_rename.h \
  src/include/dash_sched_a2av_aux.h \
  src/include/dash_sched_a2av.h \
  src/include/dash_sched_a2av_sym.h \
  src/include/dash_sched.h \
  src/include/local_generic_heap_conf.h \
  src/include/local_generic_heap.h \
  src/include/local_generic_heap_rename.h \
  src/include/prx_conf.h \
  src/include/prx.h \
  src/include/prx_rename.h \
  src/include/sl_adds.h \
  src/include/sl_adds_intern.h \
  src/include/sl_common.h \
  src/include/sl_config_global.h \
  src/include/sl_config_intern.h \
  src/include/sl_context_init.h \
  src/include/sl_context_struct.h \
  src/include/sl_data.h \
  src/include/sl_data_singles.h \
  src/include/sl_deprecated.h \
  src/include/sl_elements.h \
  src/include/sl_environment_intern.h \
  src/include/sl_globals.h \
  src/include/sl_index.h \
  src/include/sl_key.h \
  src/include/sl_macros.h \
  src/include/sl_pelem.h \
  src/include/sl_protos.h \
  src/include/sl_protos_intern.h \
  src/include/sl_protos_mpi.h \
  src/include/sl_rename.h \
  src/include/sl_rti_intern.h \
  src/include/sl_rti_tids.h \
  src/include/sl_tune_intern.h \
  src/include/sl_types.h \
  src/include/sl_types_intern.h \
  src/include/sl_xelem_call.h \
  src/include/spec_common.h \
  src/include/spec_conf.h \
  src/include/spec_core.h \
  src/include/spec_public_conf.h \
  src/include/spec_public.h \
  src/include/spec_rename.h \
  src/include/z_config.h \
  src/include/z_pack_conf.h \
  src/include/z_pack.h \
  src/include/z_pack_rename.h

sl_extra_SOURCE = \
  extra/z_pack/z_config.h \
  extra/z_pack/z_pack.c \
  extra/z_pack/z_pack_conf.h \
  extra/z_pack/z_pack.h \
  extra/z_pack/z_pack_rename.h \
  extra/zmpi_local/zmpi_local.c \
  extra/zmpi_local/zmpi_local_conf.h \
  extra/zmpi_local/zmpi_local.h \
  extra/zmpi_local/zmpi_local_rename.h \
  extra/zmpi_local/zmpil_simple.h \
  extra/zmpi_tools/zmpi_tools.c \
  extra/zmpi_tools/zmpi_tools_conf.h \
  extra/zmpi_tools/zmpi_tools.h \
  extra/zmpi_tools/zmpi_tools_rename.h \
  extra/zmpi_ataip/dash_aux_heap.h \
  extra/zmpi_ataip/dash_aux_static.h \
  extra/zmpi_ataip/dash.c \
  extra/zmpi_ataip/dash_common.h \
  extra/zmpi_ataip/dash_conf.h \
  extra/zmpi_ataip/dash_core.h \
  extra/zmpi_ataip/dash_exec.h \
  extra/zmpi_ataip/dash_exec_mpi.h \
  extra/zmpi_ataip/dash_rename.h \
  extra/zmpi_ataip/dash_sched_a2a_sym.h \
  extra/zmpi_ataip/dash_sched_a2av_aux.h \
  extra/zmpi_ataip/dash_sched_a2av.h \
  extra/zmpi_ataip/dash_sched_a2av_sym.h \
  extra/zmpi_ataip/dash_sched.h \
  extra/zmpi_ataip/local_generic_heap_conf.h \
  extra/zmpi_ataip/local_generic_heap.h \
  extra/zmpi_ataip/local_generic_heap_rename.h \
  extra/zmpi_ataip/prx_conf.h \
  extra/zmpi_ataip/prx.h \
  extra/zmpi_ataip/prx_rename.h \
  extra/zmpi_ataip/zmpi_ataip.c \
  extra/zmpi_ataip/zmpi_ataip.h \
  extra/zmpi_ataip/zmpi_ataip_rename.h \
  extra/zmpi_atasp/spec.c \
  extra/zmpi_atasp/spec_common.h \
  extra/zmpi_atasp/spec_conf.h \
  extra/zmpi_atasp/spec_core.h \
  extra/zmpi_atasp/spec_public_conf.h \
  extra/zmpi_atasp/spec_public.h \
  extra/zmpi_atasp/spec_rename.h \
  extra/zmpi_atasp/zmpi_atasp.c \
  extra/zmpi_atasp/zmpi_atasp.h \
  extra/zmpi_atasp/zmpi_atasp_rename.h

sl_extra_include_SOURCE = \
  extra/include/zmpi_ataip.h \
  extra/include/zmpi_ataip_rename.h \
  extra/include/zmpi_atasp.h \
  extra/include/zmpi_atasp_rename.h \
  extra/include/zmpi_local.h \
  extra/include/zmpi_local_rename.h \
  extra/include/zmpi_tools.h \
  extra/include/zmpi_tools_rename.h

libfcs_pepc_pepckeys_la_SOURCES = \
  config_pepc_sort.h \
  fortran2c_types.h \
  mpi_common.h \
  rename_pepc_sort.h \
  sl_pepckeys/sl_config.h \
  sl_pepckeys/sl_environment.h \
  sl_pepckeys/sl_extra.h \
  sl_pepckeys/sl_tune_auto.h \
  sl_pepckeys/sl_tune.h \
  $(sl_SOURCE) \
  $(sl_extra_include_SOURCE)

libfcs_pepc_extra_la_SOURCES = $(sl_extra_SOURCE)

libfcs_pepc_mpiwrap_la_SOURCES = \
  config_pepc_sort.h \
  fortran2c_types.h \
  include/sl_pepckeys.h \
  include/sl_pepcparts.h \
  mpi_common.c \
  mpi_common.h \
  mpi_pepckeys.c \
  mpi_pepcparts.c \
  rename_pepc_sort.h \
  $(sl_extra_include_SOURCE)
