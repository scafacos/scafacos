#
#  Makefile for pepc
#

# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

AM_FCFLAGS = $(FC_MODINC)../../..
AM_CPPFLAGS = -DSVNREVISION="\"$(SVNREVISION)\"" -DSL_PEPC_PREFIX=fcs_

# -------------------------------------------------------------------------------------------

# Source files

#-------------------------------------------------------------------------------------------
# Pepc library

noinst_LIBRARIES = liblpepc.a

pepc_raw_c_sources = bgp_core.c \
                     pepc_c_utils.c \
                     pthreads_c.c \
                     atomic_ops.c
pepc_raw_f90_sources = module_atomic_ops.f90 \
                       module_base64.f90 \
                       module_box.f90 \
                       module_comm_env.f90 \
                       module_coulomb_kernels.f90 \
                       module_debug.f90 \
                       module_domains.f90 \
                       module_fmm_framework.f90 \
                       module_interaction_specific.f90 \
                       module_interaction_specific_types_XYZQ.f90 \
                       module_libpepc_main.f90 \
                       module_math_tools.f90 \
                       module_mirror_boxes.f90 \
                       module_pepc.f90 \
                       module_pepc_types.f90 \
                       module_pepc_wrappers.f90 \
                       module_sort.f90 \
                       module_spacefilling.f90 \
                       module_timings.f90 \
                       module_tree.f90 \
                       module_tree_communicator.f90 \
                       module_tree_grow.f90 \
                       module_tree_node.f90 \
                       module_utils.f90 \
                       module_vtk.f90 \
                       module_walk_pthreads.f90 \
                       pepc_scafacos.f90 \
                       pthreads_f.f90 \
                       treevars.f90

liblpepc_a_SOURCES = $(pepc_raw_c_sources)

liblpepc_a_CPPFLAGS = -I$(top_srcdir)/openpa/src

ALL_INCLUDE_HEADER = $(srcdir)/pepc_debug.h 
PP_FILES = $(pepc_raw_f90_sources)

PP_SUFFIX = pp

PP_FILES_PP.f90 = $(addsuffix .$(PP_SUFFIX).f90,$(basename $(PP_FILES))) 
PP_FILES_PP.o = $(addsuffix .$(PP_SUFFIX).o,$(basename $(PP_FILES))) 

%.$(PP_SUFFIX).f90: %.f90 $(ALL_INCLUDE_HEADER) 
	@echo "Building file '$@'" 
	cat $(ALL_INCLUDE_HEADER) $< | \
	cpp -P -C $(AM_CPPFLAGS) $(CPPFLAGS) $(DEFAULT_INCLUDES) | \
	sed '/\/\*\(.*\)\*\//{s,,!\1,p;d;};/\/\*/,/\*\//{s,/\*,,; s,\*/,,; s,^,!,p; d;}' > $@

%.o: %.f90
	@echo "target: $@, prerequisites: $^"
	$(PPFCCOMPILE) -c -o $@ $(FCFLAGS_f90) $<

liblpepc_a_LIBADD = $(PP_FILES_PP.o)

EXTRA_DIST = $(ALL_INCLUDE_HEADER) $(PP_FILES) 
CLEANFILES = $(PP_FILES_PP.f90) 

lpepc:  liblpepc.a

mods:   all

debug:	$(OBJECTS)
	@echo "Loading object files using flags $(DBFLAGS)"
	$(FC) $(LDFLAGS)  $(VISITLIBS) $(DBFLAGS)  -o $@ $(OBJECTS)

all-local: make-opa-config

make-opa-config: $(top_builddir)/openpa/src/opa_config.h
	@if [ ! -f opa_config.h ] ; then ln -s $^ opa_config.h ; fi

clean-local: clean-opa-config clean-other-files cleanps

clean-opa-config:
	rm -f opa_config.h

clean-other-files:
	$(RM) *% *~ bin/key rand *.pif core module.deps

cleanps: 
	$(RM) *.ps 

# Use generated module dependency rules.
$(eval include $(builddir)/module.deps)

tree_stats.pp.o: module_walk_communicator.pp.o

$(builddir)/module.deps: $(pepc_raw_f90_sources) $(ALL_INCLUDE_HEADER)
	SRCDIR="$(srcdir)" PRE_INCLUDE="$(ALL_INCLUDE_HEADER)" CPP="cpp $(AM_CPPFLAGS) $(CPPFLAGS) $(DEFAULT_INCLUDES)" OBJEXT="$(PP_SUFFIX).$(OBJEXT)" MODEXT="$(FC_MODEXT)" \
	$(SHELL) $(top_srcdir)/build-aux/fortran-depcomp $^ > $@
#	@echo "uname -a: `uname -a`"
#	$(MAKE) -v
#	@echo "--> $@" ; cat $@ ; echo "<-- $@"

DISTCLEANFILES = module.deps

.SECONDARY:

# prevent parallel build on nana (workaround for weird jenkins bug)
$(if $(subst nana,,$(shell uname -n)),,$(eval .NOTPARALLEL:))

#$(eval include $(srcdir)/Makefile.debug)
