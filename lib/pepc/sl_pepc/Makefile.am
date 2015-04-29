
lib_LTLIBRARIES =
noinst_LTLIBRARIES =

noinst_LTLIBRARIES += libsl_pepc.la

sl_sub_libs = libfcs_pepc_pepckeys.la libfcs_pepc_extra.la libfcs_pepc_mpiwrap.la

noinst_LTLIBRARIES += $(sl_sub_libs)

libfcs_pepc_pepckeys_la_CPPFLAGS = \
  -DSL_PREFIX=fcs_pepckeys_ -DSL_USE_MPI -I$(srcdir)/../../common/sl/include -I$(srcdir)/sl_pepckeys -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_pepc_ -I$(srcdir)/../../common/sl/extra/include

libfcs_pepc_extra_la_CPPFLAGS = \
  -DZ_PREFIX=fcs_pepc_ -DZMPI_PREFIX=fcs_pepc_ -DHAVE_Z_PACK_H -I$(srcdir)/../../common/sl/extra/z_pack -DHAVE_ZMPI_LOCAL_H -I$(srcdir)/../../common/sl/extra/zmpi_local -DHAVE_ZMPI_TOOLS_H -I$(srcdir)/../../common/sl/extra/zmpi_tools

libfcs_pepc_mpiwrap_la_CPPFLAGS = \
  -DSL_USE_MPI -I$(srcdir)/include -DNOT_sl_pepcparts -DZMPI_PREFIX=fcs_pepc_ -DHAVE_ZMPI_LOCAL_H -DHAVE_ZMPI_TOOLS_H -DHAVE_ZMPI_ATAIP_H -DHAVE_ZMPI_ATASP_H -I$(srcdir)/../../common/sl/extra/include

libsl_pepc_la_SOURCES =

libsl_pepc_la_LIBADD = $(sl_sub_libs)

sl_SOURCE = \
  ../../common/sl/base/base.c \
  ../../common/sl/base_mpi/base_mpi.c \
  ../../common/sl/base_mpi/dash.c \
  ../../common/sl/base_mpi/spec.c \
  ../../common/sl/base/z_crc32.c \
  ../../common/sl/base/z_pack.c \
  ../../common/sl/include/dash_aux_heap.h \
  ../../common/sl/include/dash_aux_static.h \
  ../../common/sl/include/dash_common.h \
  ../../common/sl/include/dash_conf.h \
  ../../common/sl/include/dash_core.h \
  ../../common/sl/include/dash_exec.h \
  ../../common/sl/include/dash_exec_sl.h \
  ../../common/sl/include/dash_rename.h \
  ../../common/sl/include/dash_sched_a2av_aux.h \
  ../../common/sl/include/dash_sched_a2av.h \
  ../../common/sl/include/dash_sched_a2av_sym.h \
  ../../common/sl/include/dash_sched.h \
  ../../common/sl/include/local_generic_heap_conf.h \
  ../../common/sl/include/local_generic_heap.h \
  ../../common/sl/include/local_generic_heap_rename.h \
  ../../common/sl/include/prx_conf.h \
  ../../common/sl/include/prx.h \
  ../../common/sl/include/prx_rename.h \
  ../../common/sl/include/sl_adds.h \
  ../../common/sl/include/sl_adds_intern.h \
  ../../common/sl/include/sl_common.h \
  ../../common/sl/include/sl_config_global.h \
  ../../common/sl/include/sl_config_intern.h \
  ../../common/sl/include/sl_context_init.h \
  ../../common/sl/include/sl_context_struct.h \
  ../../common/sl/include/sl_data.h \
  ../../common/sl/include/sl_data_singles.h \
  ../../common/sl/include/sl_deprecated.h \
  ../../common/sl/include/sl_elements.h \
  ../../common/sl/include/sl_environment_intern.h \
  ../../common/sl/include/sl_globals.h \
  ../../common/sl/include/sl_index.h \
  ../../common/sl/include/sl_key.h \
  ../../common/sl/include/sl_macros.h \
  ../../common/sl/include/sl_pelem.h \
  ../../common/sl/include/sl_protos.h \
  ../../common/sl/include/sl_protos_intern.h \
  ../../common/sl/include/sl_protos_mpi.h \
  ../../common/sl/include/sl_rename.h \
  ../../common/sl/include/sl_rti_intern.h \
  ../../common/sl/include/sl_rti_tids.h \
  ../../common/sl/include/sl_tune_intern.h \
  ../../common/sl/include/sl_types.h \
  ../../common/sl/include/sl_types_intern.h \
  ../../common/sl/include/sl_xelem_call.h \
  ../../common/sl/include/spec_common.h \
  ../../common/sl/include/spec_conf.h \
  ../../common/sl/include/spec_core.h \
  ../../common/sl/include/spec_public_conf.h \
  ../../common/sl/include/spec_public.h \
  ../../common/sl/include/spec_rename.h \
  ../../common/sl/include/z_config.h \
  ../../common/sl/include/z_pack_conf.h \
  ../../common/sl/include/z_pack.h \
  ../../common/sl/include/z_pack_rename.h

sl_extra_SOURCE = \
  ../../common/sl/extra/z_pack/z_config.h \
  ../../common/sl/extra/z_pack/z_pack.c \
  ../../common/sl/extra/z_pack/z_pack_conf.h \
  ../../common/sl/extra/z_pack/z_pack.h \
  ../../common/sl/extra/z_pack/z_pack_rename.h \
  ../../common/sl/extra/zmpi_local/zmpi_local.c \
  ../../common/sl/extra/zmpi_local/zmpi_local_conf.h \
  ../../common/sl/extra/zmpi_local/zmpi_local.h \
  ../../common/sl/extra/zmpi_local/zmpi_local_rename.h \
  ../../common/sl/extra/zmpi_local/zmpil_simple.h \
  ../../common/sl/extra/zmpi_tools/zmpi_tools.c \
  ../../common/sl/extra/zmpi_tools/zmpi_tools_conf.h \
  ../../common/sl/extra/zmpi_tools/zmpi_tools.h \
  ../../common/sl/extra/zmpi_tools/zmpi_tools_rename.h \
  ../../common/sl/extra/zmpi_ataip/dash_aux_heap.h \
  ../../common/sl/extra/zmpi_ataip/dash_aux_static.h \
  ../../common/sl/extra/zmpi_ataip/dash.c \
  ../../common/sl/extra/zmpi_ataip/dash_common.h \
  ../../common/sl/extra/zmpi_ataip/dash_conf.h \
  ../../common/sl/extra/zmpi_ataip/dash_core.h \
  ../../common/sl/extra/zmpi_ataip/dash_exec.h \
  ../../common/sl/extra/zmpi_ataip/dash_exec_mpi.h \
  ../../common/sl/extra/zmpi_ataip/dash_rename.h \
  ../../common/sl/extra/zmpi_ataip/dash_sched_a2a_sym.h \
  ../../common/sl/extra/zmpi_ataip/dash_sched_a2av_aux.h \
  ../../common/sl/extra/zmpi_ataip/dash_sched_a2av.h \
  ../../common/sl/extra/zmpi_ataip/dash_sched_a2av_sym.h \
  ../../common/sl/extra/zmpi_ataip/dash_sched.h \
  ../../common/sl/extra/zmpi_ataip/local_generic_heap_conf.h \
  ../../common/sl/extra/zmpi_ataip/local_generic_heap.h \
  ../../common/sl/extra/zmpi_ataip/local_generic_heap_rename.h \
  ../../common/sl/extra/zmpi_ataip/prx_conf.h \
  ../../common/sl/extra/zmpi_ataip/prx.h \
  ../../common/sl/extra/zmpi_ataip/prx_rename.h \
  ../../common/sl/extra/zmpi_ataip/zmpi_ataip.c \
  ../../common/sl/extra/zmpi_ataip/zmpi_ataip.h \
  ../../common/sl/extra/zmpi_ataip/zmpi_ataip_rename.h \
  ../../common/sl/extra/zmpi_atasp/spec.c \
  ../../common/sl/extra/zmpi_atasp/spec_common.h \
  ../../common/sl/extra/zmpi_atasp/spec_conf.h \
  ../../common/sl/extra/zmpi_atasp/spec_core.h \
  ../../common/sl/extra/zmpi_atasp/spec_public_conf.h \
  ../../common/sl/extra/zmpi_atasp/spec_public.h \
  ../../common/sl/extra/zmpi_atasp/spec_rename.h \
  ../../common/sl/extra/zmpi_atasp/zmpi_atasp.c \
  ../../common/sl/extra/zmpi_atasp/zmpi_atasp.h \
  ../../common/sl/extra/zmpi_atasp/zmpi_atasp_rename.h

sl_extra_include_SOURCE = \
  ../../common/sl/extra/include/zmpi_ataip.h \
  ../../common/sl/extra/include/zmpi_ataip_rename.h \
  ../../common/sl/extra/include/zmpi_atasp.h \
  ../../common/sl/extra/include/zmpi_atasp_rename.h \
  ../../common/sl/extra/include/zmpi_local.h \
  ../../common/sl/extra/include/zmpi_local_rename.h \
  ../../common/sl/extra/include/zmpi_tools.h \
  ../../common/sl/extra/include/zmpi_tools_rename.h

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
