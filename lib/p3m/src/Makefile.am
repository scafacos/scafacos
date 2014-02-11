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
SUBDIRS = . tests

noinst_LTLIBRARIES = libp3m.la
libp3m_la_CPPFLAGS = $(fftw3_CPPFLAGS) -I $(top_srcdir)/src \
	-I $(top_srcdir)/lib -I $(top_srcdir)/lib/common/fcs-common
libp3m_la_SOURCES = \
	p3mconfig.hpp \
	utils.hpp \
	Communication.hpp Communication.cpp \
	caf.cpp caf.hpp \
	Parallel3DFFT.cpp Parallel3DFFT.hpp \
	influence_function.cpp influence_function.hpp \
	init.cpp init.hpp \
	prepare.cpp prepare.hpp \
	timing.cpp timing.hpp \
	tune.cpp tune.hpp \
	tune_broadcast.cpp tune_broadcast.hpp \
	run.cpp run.hpp \
	types.hpp \
	error_estimate.cpp error_estimate.hpp
