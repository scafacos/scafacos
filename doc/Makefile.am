# Copyright (C) 2011,2012,2013 The ScaFaCoS project
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
	figures/logo.png

EXTRA_DIST = $(manual_FILES)
CLEANFILES =
MOSTLYCLEANFILES =
CLEAN_LOCAL =
doc_DATA =

manual: manual.pdf

if HAVE_LATEX
manual.pdf: $(manual_FILES) latex-manual.sh latexinc.sh
	sh latex-manual.sh

doc_DATA   += manual.pdf
EXTRA_DIST += manual.pdf
CLEANFILES += manual.pdf
MOSTLYCLEANFILES += \
	*.aux *.aux.bak \
	manual.bbl manual.blg \
	manual.idx manual.idx.bak manual.ilg manual.ind \
	manual.log manual.out manual.toc manual.lof manual.lot manual.tdo
else
manual.pdf: 
	@echo "No complete LaTeX-installation was found in your PATH."
	@echo "Can't build the User's Guide without pdflatex, makeindex and bibtex."
	@echo "Install these programs and rerun configure."
endif

html: doxyfile.stamp
manpage: doxyfile.stamp

if HAVE_DOXYGEN
doxyfile.stamp: Doxyfile
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

EXTRA_DIST += html man
CLEANFILES += doxyfile.stamp
CLEAN_LOCAL += html man
    
# I assume, that those should work but dont - that is why we add the two rules below
#html_DATA = html
#man_MANS  = man/man3/*

install-data-local:	html
	@$(NORMAL_INSTALL)
	$(MKDIR_P) $(htmldir)/html
	$(INSTALL_DATA) $(builddir)/html/* $(htmldir)/html/
	$(MKDIR_P) $(mandir)/man3
	$(INSTALL_DATA) $(builddir)/man/man3/* $(mandir)/man3/

## Warning: the following rule assumes that this package is the
## sole owner of files below $(htmldir)/html!
uninstall-local:
	rm $(htmldir)/html/*
	-rmdir $(htmldir)/html

else
doxyfile.stamp: 
	@echo "No Doxygen was found in your PATH."
	@echo "Can't build the HTML and manpage documentation without Doxygen."
	@echo "Install this program and rerun configure."
endif

doc: manual html manpage

UPLOADTMP=$(builddir)/upload.tmp
doc-upload: manual html
	git clone https://github.com/scafacos/scafacos.github.io $(UPLOADTMP)
	cd $(UPLOADTMP); cp ../manual.pdf ./files/pub/libfcs_manual.pdf; git add ./files/pub/libfcs_manual.pdf
	cd $(UPLOADTMP); git rm -r ./doxygen/*; cp -r ../html ./doxygen; git add ./doxygen/*
	-cd $(UPLOADTMP); git commit; git push origin master
	-rm -rf $(UPLOADTMP)

.PHONY: doc manual html manpage doc-upload

clean-local:
	-rm -rf $(CLEAN_LOCAL)
