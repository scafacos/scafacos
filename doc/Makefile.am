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
#


# add all source files here
manual_FILES = \
	manual.tex \
	introduction.tex \
	compilation.tex \
	interface.tex \
	fmm.tex \
	memd.tex \
	mmm1d.tex \
	mmm2d.tex \
	p2nfft.tex \
	p3m.tex \
	pepc.tex \
	pp3mg.tex \
	vmg.tex \
	direct.tex \
	ewald.tex \
	tests.tex \
	buildsystem.tex \
	implementation.tex \
	functionlist.tex \
	todonotes.sty \
	bibliography.bib \
	figures/hball_500.pdf \
	figures/hballneg_500.pdf \
	figures/p2nfft_ewald_split.pdf \
	figures/vmg_short-long-splitting.tex figures/vmg_short-long-splitting.pdf \
    figures/logo.png \
    Doxyfile

EXTRA_DIST = $(manual_FILES)
if HAVE_LATEX
EXTRA_DIST += manual.pdf
endif

.PHONY: doc manual html

doc: manual doxyfile.stamp
manual: manual.pdf
manual.pdf: $(manual_FILES) latex-manual.sh latexinc.sh
if HAVE_LATEX
	sh latex-manual.sh
else
	@echo "No complete LaTeX-installation was found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these and rerun configure."
endif
    
if HAVE_DOXYGEN
mandir = $(top_srcdir)/doc/man/man3/

dist_man_MANS = $(mandir)/man_page_1.3 $(mandir)/man_page_2.3
$(mandir)/man_page_1.3: doxyfile.stamp
$(mandir)/man_page_2.3: doxyfile.stamp

doxyfile.stamp:
   $(DOXYGEN) Doxyfile
   echo Timestamp > doxyfile.stamp
endif

CLEANFILES = manual.pdf html $(top_srcdir)/doc/man
if HAVE_DOXYGEN
CLEANFILES += doxyfile.stamp
endif

MOSTLYCLEANFILES = \
	*.aux *.aux.bak\
	manual.bbl manual.blg \
	manual.idx manual.idx.bak manual.ilg manual.ind\
	manual.log manual.out manual.toc manual.lof manual.lot manual.tdo
