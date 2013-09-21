
lib_LIBRARIES =
noinst_LIBRARIES =

if ENABLE_SINGLE_LIB
noinst_LIBRARIES += libfcs_gridsort.a
else
lib_LIBRARIES += libfcs_gridsort.a
endif

sl_sub_libs = libsl_back__p.a libsl_back_f_.a libsl_back_fp.a libsl_back_x.a libsl_forw.a libsl_extra.a libsl_mpiwrap.a

noinst_LIBRARIES += $(sl_sub_libs)

srcdir_sl=$(srcdir)
srcdir_sl_source=$(srcdir)/../sl
srcdir_sl_extra=$(srcdir)/../sl/extra

libsl_back__p_a_CPPFLAGS = \
  -DSL_PREFIX=fcs_back__p_ -DSL_USE_MPI -I$(srcdir_sl_source)/include -I$(srcdir_sl)/sl_back__p -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_gridsort_ -I$(srcdir_sl_extra)/include

libsl_back_f__a_CPPFLAGS = \
  -DSL_PREFIX=fcs_back_f__ -DSL_USE_MPI -I$(srcdir_sl_source)/include -I$(srcdir_sl)/sl_back_f_ -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_gridsort_ -I$(srcdir_sl_extra)/include

libsl_back_fp_a_CPPFLAGS = \
  -DSL_PREFIX=fcs_back_fp_ -DSL_USE_MPI -I$(srcdir_sl_source)/include -I$(srcdir_sl)/sl_back_fp -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_gridsort_ -I$(srcdir_sl_extra)/include

libsl_back_x_a_CPPFLAGS = \
  -DSL_PREFIX=fcs_back_x_ -DSL_USE_MPI -I$(srcdir_sl_source)/include -I$(srcdir_sl)/sl_back_x -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_gridsort_ -I$(srcdir_sl_extra)/include

libsl_forw_a_CPPFLAGS = \
  -DSL_PREFIX=fcs_forw_ -DSL_USE_MPI -I$(srcdir_sl_source)/include -I$(srcdir_sl)/sl_forw -DSL_EXTRA -DSL_EXTRA_PREFIX=fcs_gridsort_ -I$(srcdir_sl_extra)/include

libsl_extra_a_CPPFLAGS = \
  -DZ_PREFIX=fcs_gridsort_ -DZMPI_PREFIX=fcs_gridsort_ -DHAVE_Z_PACK_H -I$(srcdir_sl_extra)/z_pack -DHAVE_ZMPI_LOCAL_H -I$(srcdir_sl_extra)/zmpi_local -DHAVE_ZMPI_TOOLS_H -I$(srcdir_sl_extra)/zmpi_tools

libsl_mpiwrap_a_CPPFLAGS = \
  -DSL_USE_MPI -I$(srcdir_sl)/include -I$(top_srcdir)/lib -DZMPI_PREFIX=fcs_gridsort_ -DHAVE_ZMPI_LOCAL_H -DHAVE_ZMPI_TOOLS_H -DHAVE_ZMPI_ATAIP_H -DHAVE_ZMPI_ATASP_H -I$(srcdir_sl_extra)/include

libfcs_gridsort_a_SOURCES =

libfcs_gridsort.a: $(sl_sub_libs)
	rm -f $@
	$(AR) $(ARFLAGS) $@ *.o
	$(RANLIB) $@

sl_SOURCE = \
  $(srcdir_sl_source)/base/base.c \
  $(srcdir_sl_source)/base_mpi/base_mpi.c \
  $(srcdir_sl_source)/base_mpi/dash.c \
  $(srcdir_sl_source)/base_mpi/spec.c \
  $(srcdir_sl_source)/base/z_crc32.c \
  $(srcdir_sl_source)/base/z_pack.c \
  $(srcdir_sl_source)/include/dash_aux_heap.h \
  $(srcdir_sl_source)/include/dash_aux_static.h \
  $(srcdir_sl_source)/include/dash_common.h \
  $(srcdir_sl_source)/include/dash_conf.h \
  $(srcdir_sl_source)/include/dash_core.h \
  $(srcdir_sl_source)/include/dash_exec.h \
  $(srcdir_sl_source)/include/dash_exec_sl.h \
  $(srcdir_sl_source)/include/dash_rename.h \
  $(srcdir_sl_source)/include/dash_sched_a2av_aux.h \
  $(srcdir_sl_source)/include/dash_sched_a2av.h \
  $(srcdir_sl_source)/include/dash_sched_a2av_sym.h \
  $(srcdir_sl_source)/include/dash_sched.h \
  $(srcdir_sl_source)/include/local_generic_heap_conf.h \
  $(srcdir_sl_source)/include/local_generic_heap.h \
  $(srcdir_sl_source)/include/local_generic_heap_rename.h \
  $(srcdir_sl_source)/include/prx_conf.h \
  $(srcdir_sl_source)/include/prx.h \
  $(srcdir_sl_source)/include/prx_rename.h \
  $(srcdir_sl_source)/include/sl_adds.h \
  $(srcdir_sl_source)/include/sl_adds_intern.h \
  $(srcdir_sl_source)/include/sl_common.h \
  $(srcdir_sl_source)/include/sl_config_global.h \
  $(srcdir_sl_source)/include/sl_config_intern.h \
  $(srcdir_sl_source)/include/sl_context_init.h \
  $(srcdir_sl_source)/include/sl_context_struct.h \
  $(srcdir_sl_source)/include/sl_data.h \
  $(srcdir_sl_source)/include/sl_data_singles.h \
  $(srcdir_sl_source)/include/sl_deprecated.h \
  $(srcdir_sl_source)/include/sl_elements.h \
  $(srcdir_sl_source)/include/sl_environment_intern.h \
  $(srcdir_sl_source)/include/sl_globals.h \
  $(srcdir_sl_source)/include/sl_index.h \
  $(srcdir_sl_source)/include/sl_key.h \
  $(srcdir_sl_source)/include/sl_pelem.h \
  $(srcdir_sl_source)/include/sl_protos.h \
  $(srcdir_sl_source)/include/sl_protos_mpi.h \
  $(srcdir_sl_source)/include/sl_rename.h \
  $(srcdir_sl_source)/include/sl_rti_intern.h \
  $(srcdir_sl_source)/include/sl_rti_tids.h \
  $(srcdir_sl_source)/include/sl_tune_intern.h \
  $(srcdir_sl_source)/include/sl_types.h \
  $(srcdir_sl_source)/include/sl_types_intern.h \
  $(srcdir_sl_source)/include/sl_xelem_call.h \
  $(srcdir_sl_source)/include/spec_conf.h \
  $(srcdir_sl_source)/include/spec_core.h \
  $(srcdir_sl_source)/include/spec_public_conf.h \
  $(srcdir_sl_source)/include/spec_public.h \
  $(srcdir_sl_source)/include/spec_rename.h \
  $(srcdir_sl_source)/include/z_config.h \
  $(srcdir_sl_source)/include/z_pack_conf.h \
  $(srcdir_sl_source)/include/z_pack.h \
  $(srcdir_sl_source)/include/z_pack_rename.h \
  $(srcdir_sl_extra)/include/zmpi_ataip.h \
  $(srcdir_sl_extra)/include/zmpi_ataip_rename.h \
  $(srcdir_sl_extra)/include/zmpi_atasp.h \
  $(srcdir_sl_extra)/include/zmpi_atasp_rename.h \
  $(srcdir_sl_extra)/include/zmpi_local.h \
  $(srcdir_sl_extra)/include/zmpi_local_rename.h \
  $(srcdir_sl_extra)/include/zmpi_tools.h \
  $(srcdir_sl_extra)/include/zmpi_tools_rename.h

libsl_back__p_a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/z_tools.h \
  $(srcdir_sl)/sl_back__p/sl_config.h \
  $(srcdir_sl)/sl_back__p/sl_environment.h \
  $(srcdir_sl)/sl_back__p/sl_extra.h \
  $(srcdir_sl)/sl_back__p/sl_tune_auto.h \
  $(srcdir_sl)/sl_back__p/sl_tune.h \
  $(sl_SOURCE)

libsl_back_f__a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/z_tools.h \
  $(srcdir_sl)/sl_back_f_/sl_config.h \
  $(srcdir_sl)/sl_back_f_/sl_environment.h \
  $(srcdir_sl)/sl_back_f_/sl_extra.h \
  $(srcdir_sl)/sl_back_f_/sl_tune_auto.h \
  $(srcdir_sl)/sl_back_f_/sl_tune.h \
  $(sl_SOURCE)

libsl_back_fp_a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/z_tools.h \
  $(srcdir_sl)/sl_back_fp/sl_config.h \
  $(srcdir_sl)/sl_back_fp/sl_environment.h \
  $(srcdir_sl)/sl_back_fp/sl_extra.h \
  $(srcdir_sl)/sl_back_fp/sl_tune_auto.h \
  $(srcdir_sl)/sl_back_fp/sl_tune.h \
  $(sl_SOURCE)

libsl_back_x_a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/z_tools.h \
  $(srcdir_sl)/sl_back_x/sl_config.h \
  $(srcdir_sl)/sl_back_x/sl_environment.h \
  $(srcdir_sl)/sl_back_x/sl_extra.h \
  $(srcdir_sl)/sl_back_x/sl_tune_auto.h \
  $(srcdir_sl)/sl_back_x/sl_tune.h \
  $(sl_SOURCE)

libsl_forw_a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/z_tools.h \
  $(srcdir_sl)/sl_forw/sl_config.h \
  $(srcdir_sl)/sl_forw/sl_environment.h \
  $(srcdir_sl)/sl_forw/sl_extra.h \
  $(srcdir_sl)/sl_forw/sl_tune_auto.h \
  $(srcdir_sl)/sl_forw/sl_tune.h \
  $(sl_SOURCE)

libsl_extra_a_SOURCES = \
  $(srcdir_sl_extra)/z_pack/z_config.h \
  $(srcdir_sl_extra)/z_pack/z_pack.c \
  $(srcdir_sl_extra)/z_pack/z_pack_conf.h \
  $(srcdir_sl_extra)/z_pack/z_pack.h \
  $(srcdir_sl_extra)/z_pack/z_pack_rename.h \
  $(srcdir_sl_extra)/zmpi_local/zmpi_local.c \
  $(srcdir_sl_extra)/zmpi_local/zmpi_local_conf.h \
  $(srcdir_sl_extra)/zmpi_local/zmpi_local.h \
  $(srcdir_sl_extra)/zmpi_local/zmpi_local_rename.h \
  $(srcdir_sl_extra)/zmpi_local/zmpil_simple.h \
  $(srcdir_sl_extra)/zmpi_tools/zmpi_tools.c \
  $(srcdir_sl_extra)/zmpi_tools/zmpi_tools_conf.h \
  $(srcdir_sl_extra)/zmpi_tools/zmpi_tools.h \
  $(srcdir_sl_extra)/zmpi_tools/zmpi_tools_rename.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_aux_heap.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_aux_static.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash.c \
  $(srcdir_sl_extra)/zmpi_ataip/dash_common.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_conf.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_core.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_exec.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_exec_mpi.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_rename.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_sched_a2a_sym.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_sched_a2av_aux.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_sched_a2av.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_sched_a2av_sym.h \
  $(srcdir_sl_extra)/zmpi_ataip/dash_sched.h \
  $(srcdir_sl_extra)/zmpi_ataip/local_generic_heap_conf.h \
  $(srcdir_sl_extra)/zmpi_ataip/local_generic_heap.h \
  $(srcdir_sl_extra)/zmpi_ataip/local_generic_heap_rename.h \
  $(srcdir_sl_extra)/zmpi_ataip/prx_conf.h \
  $(srcdir_sl_extra)/zmpi_ataip/prx.h \
  $(srcdir_sl_extra)/zmpi_ataip/prx_rename.h \
  $(srcdir_sl_extra)/zmpi_ataip/zmpi_ataip.c \
  $(srcdir_sl_extra)/zmpi_ataip/zmpi_ataip.h \
  $(srcdir_sl_extra)/zmpi_ataip/zmpi_ataip_rename.h \
  $(srcdir_sl_extra)/zmpi_atasp/spec.c \
  $(srcdir_sl_extra)/zmpi_atasp/spec_conf.h \
  $(srcdir_sl_extra)/zmpi_atasp/spec_core.h \
  $(srcdir_sl_extra)/zmpi_atasp/spec_public_conf.h \
  $(srcdir_sl_extra)/zmpi_atasp/spec_public.h \
  $(srcdir_sl_extra)/zmpi_atasp/spec_rename.h \
  $(srcdir_sl_extra)/zmpi_atasp/zmpi_atasp.c \
  $(srcdir_sl_extra)/zmpi_atasp/zmpi_atasp.h \
  $(srcdir_sl_extra)/zmpi_atasp/zmpi_atasp_rename.h

libsl_mpiwrap_a_SOURCES = \
  $(srcdir_sl)/bounds_xyz2coords.h \
  $(srcdir_sl)/common.h \
  $(srcdir_sl)/gridsort.c \
  $(srcdir_sl)/gridsort_front_tproc.h \
  $(srcdir_sl)/gridsort.h \
  $(srcdir_sl)/gridsort_resort.c \
  $(srcdir_sl)/gridsort_resort.h \
  $(srcdir_sl)/include/sl_back_f_.h \
  $(srcdir_sl)/include/sl_back_fp.h \
  $(srcdir_sl)/include/sl_back__p.h \
  $(srcdir_sl)/include/sl_back_x.h \
  $(srcdir_sl)/include/sl_forw.h \
  $(srcdir_sl)/z_tools.h
