EXTRA_DIST = versamg.doxygen

html-local:
	cd $(srcdir) && $(DOXYGEN) versamg.doxygen

clean-local: clean-html
clean-html:
	-cd $(srcdir) && rm -rf html latex
