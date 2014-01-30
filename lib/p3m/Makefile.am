# Copyright (C) 2012,2013,2014 Olaf Lenz
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
SUBDIRS = src .

if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_p3m.la
else
lib_LTLIBRARIES = libfcs_p3m.la
endif

libfcs_p3m_la_CPPFLAGS = -I$(top_srcdir)/src $(fftw3_CPPFLAGS)
libfcs_p3m_la_SOURCES = \
	scafacos.h \
	scafacos.cpp run.cpp
libfcs_p3m_la_LIBADD = src/libp3m.la

