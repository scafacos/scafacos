# Copyright (C) 2014 Olaf Lenz
#
# This file is part of ScaFaCoS.
#
# ScaFaCoS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ScaFaCoS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
DUMMY:=$(eval include $(top_builddir)/package/scafacos-fcs.mk)

check_PROGRAMS = compute_error
compute_error_SOURCES = compute_error.cpp
compute_error_LDADD = $(SCAFACOS_MK_LDADD) $(top_builddir)/src/libfcs.la
compute_error_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
compute_error_CPPFLAGS = -I $(top_srcdir)/src -I $(srcdir)/.. $(fftw3_CPPFLAGS)
