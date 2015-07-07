PYTHONFILES = setup.py scafacos.pyx test.py

EXTRA_DIST = $(COPYFILES)

# install python files in a VPATH build.
all-local: make-python-files
make-python-files:
	-test . != "$(srcdir)" && for file in $(PYTHONFILES); do \
	  $(INSTALL_DATA) "$(srcdir)/$$file" .; \
	done

distclean-local: clean-python-files
clean-python-files:
	-test . != "$(srcdir)" && rm -f $(PYTHONFILES)

.PHONY: make-python-files clean-python-files
