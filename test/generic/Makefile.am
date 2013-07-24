# Copyright (C) 2011 The ScaFaCoS project
#  
# This file is part of ScaFaCoS.
#  
# ScaFaCoS is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
#  ScaFaCoS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 

# Get linker flags from scafacos-fcs.mk. (The assignment and evaluation is
# required to force Automake to include 'scafacos-fcs.mk' early enough.)
DUMMY:=$(eval include $(top_builddir)/package/scafacos-fcs.mk)

# Preprocessor flags.
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib

# Libraries.
LDADD = $(SCAFACOS_MK_LDADD)

if HAVE_IOMANIP
check_PROGRAMS = scafacos_test
scafacos_test_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
scafacos_test_SOURCES = \
	scafacos_test.cpp \
	Testcase.cpp Testcase.hpp \
	Generator.cpp Generator.hpp \
	Integration.cpp Integration.hpp \
	PointSequence.hpp \
	common.cpp common.hpp \
	rapidxml/rapidxml.hpp \
	rapidxml/rapidxml_iterators.hpp \
	rapidxml/rapidxml_print.hpp \
	rapidxml/rapidxml_utils.hpp
endif

EXTRA_DIST = \
	systems/xyz2xml.py \
	systems/xml2xyz.py \
	systems/xml2vtf.py \
	systems/xml2vtk.py \
	scafacos_test_all.sh \
	tune_alpha_2dp.py \
	tune_alpha_3dp.py \
	tune_alpha_cloud_wall_3dp.py \
	opt_alpha_cloud_wall_2dp.py \
	opt_alpha_cloud_wall_3dp.py \
	compute_nonperiodic_references.sh \
	rapidxml/manual.html \
	rapidxml/license.txt

EXTRA_DIST += \
	systems/3d-periodic/cloud_wall_8100.xml.gz \
	systems/3d-periodic/cloud_wall_102900.xml.gz \
	systems/3d-periodic/cloud_wall_1012500.xml.gz \
	systems/3d-periodic/silica_melt_12960.xml.gz \
	systems/3d-periodic/silica_melt_103680.xml.gz \
	systems/3d-periodic/silica_melt_829440.xml.gz \
	systems/nonperiodic/cloud_wall_8100.xml.gz \
	systems/nonperiodic/silica_melt_12960.xml.gz

# Generate symbolic link to systems subdir
all-local: make-data-symlink
make-data-symlink:
	-@if test "$(srcdir)" != . && test ! -r systems; then \
	  ln -s "$(srcdir)/systems" .; \
	fi
	-@if test "$(srcdir)" != . && test ! -r hammersley_generator; then \
	  ln -s "$(srcdir)/hammersley_generator" .; \
	fi

distclean-local: clean-data-symlink
clean-data-symlink:
	-@if test "$(srcdir)" != .; then \
	  rm -f systems hammersley_generator; \
	fi

.PHONY: make-data-symlink clean-data-symlink

# Distribute all input data, but not any SVN files.
#dist-hook:
#	cp -pR $(srcdir)/systems $(distdir)/systems
#	rm -rf `find $(distdir)/systems -name .svn`
