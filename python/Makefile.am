COPYFILES = setup.py scafacos.pyx test.py

EXTRA_DIST = $(COPYFILES)

.PHONY: copy-files
FORCE:
copy-files: FORCE
	for file in $(COPYFILES); do \
	  if test ! -e $$file || test ! $$file -ef $(abs_srcdir)/$$file; then \
	    $(INSTALL_DATA) $(abs_srcdir)/$$file $(abs_builddir); \
	  fi; \
	done

all-local: copy-files

