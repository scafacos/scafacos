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

# Get linker flags from scafacos.mk generated during configure.
# * The 'indirection' of the @...@ substitution allows to inject an include
#   command into the Makefile without being expanded by Automake.
# * The variable definition ensures that Automake adds this line before the
#   build rules.
# * Using the command as a variable name 'activates' the deferred expansion of
#   the SCAFACOS_MK_LDADD and SCAFACOS_MK_DEPS variables from scafacos.mk.
DUMMY_SCAFACOS_MK_@EVAL_INCLUDE_PACKAGE_SCAFACOS_MK@ = DUMMY

check_PROGRAMS = compute_error
compute_error_SOURCES = compute_error.cpp
compute_error_LDADD = ../libp3m.la 
compute_error_CPPFLAGS = -I $(top_srcdir)/src -I $(srcdir)/.. $(fftw3_CPPFLAGS)
