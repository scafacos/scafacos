
lib_LTLIBRARIES =
noinst_LTLIBRARIES =

if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES += libfcs_resort.la
else
lib_LTLIBRARIES += libfcs_resort.la
endif

sl_sub_libs = libfcs_resort_extra.la libfcs_resort_mpiwrap.la

noinst_LTLIBRARIES += $(sl_sub_libs)

libfcs_resort_extra_la_CPPFLAGS = \
  -DZ_PREFIX=fcs_resort_ -DZMPI_PREFIX=fcs_resort_ -DHAVE_Z_PACK_H -I$(srcdir)/../sl/extra/z_pack -DHAVE_ZMPI_LOCAL_H -I$(srcdir)/../sl/extra/zmpi_local -DHAVE_ZMPI_TOOLS_H -I$(srcdir)/../sl/extra/zmpi_tools

libfcs_resort_mpiwrap_la_CPPFLAGS = \
  -DSL_USE_MPI -I$(top_srcdir)/lib -I$(srcdir)/include -DZMPI_PREFIX=fcs_resort_ -DHAVE_ZMPI_LOCAL_H -DHAVE_ZMPI_TOOLS_H -DHAVE_ZMPI_ATAIP_H -DHAVE_ZMPI_ATASP_H -I$(srcdir)/../sl/extra/include

libfcs_resort_la_SOURCES =

libfcs_resort_la_LIBADD = $(sl_sub_libs)

sl_SOURCE = \
  ../sl/base/base.c \
  ../sl/base_mpi/base_mpi.c \
  ../sl/base_mpi/dash.c \
  ../sl/base_mpi/spec.c \
  ../sl/base/z_crc32.c \
  ../sl/base/z_pack.c \
  ../sl/include/dash_aux_heap.h \
  ../sl/include/dash_aux_static.h \
  ../sl/include/dash_common.h \
  ../sl/include/dash_conf.h \
  ../sl/include/dash_core.h \
  ../sl/include/dash_exec.h \
  ../sl/include/dash_exec_sl.h \
  ../sl/include/dash_rename.h \
  ../sl/include/dash_sched_a2av_aux.h \
  ../sl/include/dash_sched_a2av.h \
  ../sl/include/dash_sched_a2av_sym.h \
  ../sl/include/dash_sched.h \
  ../sl/include/local_generic_heap_conf.h \
  ../sl/include/local_generic_heap.h \
  ../sl/include/local_generic_heap_rename.h \
  ../sl/include/prx_conf.h \
  ../sl/include/prx.h \
  ../sl/include/prx_rename.h \
  ../sl/include/sl_adds.h \
  ../sl/include/sl_adds_intern.h \
  ../sl/include/sl_common.h \
  ../sl/include/sl_config_global.h \
  ../sl/include/sl_config_intern.h \
  ../sl/include/sl_context_init.h \
  ../sl/include/sl_context_struct.h \
  ../sl/include/sl_data.h \
  ../sl/include/sl_data_singles.h \
  ../sl/include/sl_deprecated.h \
  ../sl/include/sl_elements.h \
  ../sl/include/sl_environment_intern.h \
  ../sl/include/sl_globals.h \
  ../sl/include/sl_index.h \
  ../sl/include/sl_key.h \
  ../sl/include/sl_macros.h \
  ../sl/include/sl_pelem.h \
  ../sl/include/sl_protos.h \
  ../sl/include/sl_protos_intern.h \
  ../sl/include/sl_protos_mpi.h \
  ../sl/include/sl_rename.h \
  ../sl/include/sl_rti_intern.h \
  ../sl/include/sl_rti_tids.h \
  ../sl/include/sl_tune_intern.h \
  ../sl/include/sl_types.h \
  ../sl/include/sl_types_intern.h \
  ../sl/include/sl_xelem_call.h \
  ../sl/include/spec_common.h \
  ../sl/include/spec_conf.h \
  ../sl/include/spec_core.h \
  ../sl/include/spec_public_conf.h \
  ../sl/include/spec_public.h \
  ../sl/include/spec_rename.h \
  ../sl/include/z_config.h \
  ../sl/include/z_pack_conf.h \
  ../sl/include/z_pack.h \
  ../sl/include/z_pack_rename.h

sl_extra_SOURCE = \
  ../sl/extra/z_pack/z_config.h \
  ../sl/extra/z_pack/z_pack.c \
  ../sl/extra/z_pack/z_pack_conf.h \
  ../sl/extra/z_pack/z_pack.h \
  ../sl/extra/z_pack/z_pack_rename.h \
  ../sl/extra/zmpi_local/zmpi_local.c \
  ../sl/extra/zmpi_local/zmpi_local_conf.h \
  ../sl/extra/zmpi_local/zmpi_local.h \
  ../sl/extra/zmpi_local/zmpi_local_rename.h \
  ../sl/extra/zmpi_local/zmpil_simple.h \
  ../sl/extra/zmpi_tools/zmpi_tools.c \
  ../sl/extra/zmpi_tools/zmpi_tools_conf.h \
  ../sl/extra/zmpi_tools/zmpi_tools.h \
  ../sl/extra/zmpi_tools/zmpi_tools_rename.h \
  ../sl/extra/zmpi_ataip/dash_aux_heap.h \
  ../sl/extra/zmpi_ataip/dash_aux_static.h \
  ../sl/extra/zmpi_ataip/dash.c \
  ../sl/extra/zmpi_ataip/dash_common.h \
  ../sl/extra/zmpi_ataip/dash_conf.h \
  ../sl/extra/zmpi_ataip/dash_core.h \
  ../sl/extra/zmpi_ataip/dash_exec.h \
  ../sl/extra/zmpi_ataip/dash_exec_mpi.h \
  ../sl/extra/zmpi_ataip/dash_rename.h \
  ../sl/extra/zmpi_ataip/dash_sched_a2a_sym.h \
  ../sl/extra/zmpi_ataip/dash_sched_a2av_aux.h \
  ../sl/extra/zmpi_ataip/dash_sched_a2av.h \
  ../sl/extra/zmpi_ataip/dash_sched_a2av_sym.h \
  ../sl/extra/zmpi_ataip/dash_sched.h \
  ../sl/extra/zmpi_ataip/local_generic_heap_conf.h \
  ../sl/extra/zmpi_ataip/local_generic_heap.h \
  ../sl/extra/zmpi_ataip/local_generic_heap_rename.h \
  ../sl/extra/zmpi_ataip/prx_conf.h \
  ../sl/extra/zmpi_ataip/prx.h \
  ../sl/extra/zmpi_ataip/prx_rename.h \
  ../sl/extra/zmpi_ataip/zmpi_ataip.c \
  ../sl/extra/zmpi_ataip/zmpi_ataip.h \
  ../sl/extra/zmpi_ataip/zmpi_ataip_rename.h \
  ../sl/extra/zmpi_atasp/spec.c \
  ../sl/extra/zmpi_atasp/spec_common.h \
  ../sl/extra/zmpi_atasp/spec_conf.h \
  ../sl/extra/zmpi_atasp/spec_core.h \
  ../sl/extra/zmpi_atasp/spec_public_conf.h \
  ../sl/extra/zmpi_atasp/spec_public.h \
  ../sl/extra/zmpi_atasp/spec_rename.h \
  ../sl/extra/zmpi_atasp/zmpi_atasp.c \
  ../sl/extra/zmpi_atasp/zmpi_atasp.h \
  ../sl/extra/zmpi_atasp/zmpi_atasp_rename.h

sl_extra_include_SOURCE = \
  ../sl/extra/include/zmpi_ataip.h \
  ../sl/extra/include/zmpi_ataip_rename.h \
  ../sl/extra/include/zmpi_atasp.h \
  ../sl/extra/include/zmpi_atasp_rename.h \
  ../sl/extra/include/zmpi_local.h \
  ../sl/extra/include/zmpi_local_rename.h \
  ../sl/extra/include/zmpi_tools.h \
  ../sl/extra/include/zmpi_tools_rename.h

libfcs_resort_extra_la_SOURCES = $(sl_extra_SOURCE)

libfcs_resort_mpiwrap_la_SOURCES = \
  resort.c \
  resort.h \
  z_tools.h \
  $(sl_extra_include_SOURCE)
