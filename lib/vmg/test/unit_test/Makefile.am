if COND_HAVE_BOOST_UNIT_TEST

check_PROGRAMS += \
	boost_test

boost_test_SOURCES = \
	unit_test/boost_test.cpp \
	unit_test/library/dirichlet_cs.cpp \
	unit_test/library/dirichlet_cs_mpi.cpp \
	unit_test/library/dirichlet_fas.cpp \
	unit_test/library/dirichlet_fas_lr.cpp \
	unit_test/library/dirichlet_fas_lr_mpi.cpp \
	unit_test/library/dirichlet_fas_mpi.cpp \
	unit_test/library/force_sinus.cpp \
	unit_test/library/periodic_cs.cpp \
	unit_test/library/periodic_cs_mpi.cpp \
	unit_test/library/periodic_fas.cpp \
	unit_test/library/periodic_fas_mpi.cpp \
	unit_test/unit_test/interpolate_test.cpp \
	unit_test/unit_test/polynomial_test.cpp \
	unit_test/unit_test/smoother_test.cpp \
	unit_test/unit_test/solver_test.cpp

boost_test_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/boost_unit_test -I$(top_srcdir)/test/interfaces $(VTK_CXXFLAGS) $(BOOST_CPPFLAGS)
boost_test_LDADD = \
	$(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
	$(LIBVMG) \
	libinterfaces.la \
	$(VTK_LDFLAGS) \
	$(VTK_LIBS) \
	$(BOOST_LDFLAGS)

if HAVE_LAPACK
boost_test_LDADD += \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) \
	$(FCLIBS)
endif

endif
