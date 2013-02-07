#! /bin/sh
############################################################################
#
# NOTE: If you just want to build FFTW, do not use this file.  Just use
# the ordinary ./configure && make commmands as described in the installation
# section of the manual.
#
# This file is only for users that want to generate their own codelets,
# as described in the "generating your own code" section of the manual.
#
############################################################################

# touch ChangeLog
# 
# echo "PLEASE IGNORE WARNINGS AND ERRORS"
# 
# # paranoia: sometimes autoconf doesn't get things right the first time
# rm -rf autom4te.cache
# autoreconf --verbose --install --symlink --force
# autoreconf --verbose --install --symlink --force
# autoreconf --verbose --install --symlink --force
# 
# rm -f config.cache
# 
# # --enable-maintainer-mode enables build of genfft and automatic
# # rebuild of codelets whenever genfft changes
# (
#     ./configure --disable-shared --enable-maintainer-mode --enable-threads $*
# )

# end of paranoia:

# alias to allow for systems having glibtoolize
#alias libtoolize=$(type -p glibtoolize libtoolize | head -1)
alias libtoolize=$(for l in glibtoolize libtoolize ; do which $l ; done | head -1)

touch ChangeLog

echo "PLEASE IGNORE WARNINGS AND ERRORS"

rm -rf autom4te.cache
libtoolize
autoreconf --verbose --install --force

rm -f config.cache

# Add dependency tracking support for IBM C/C++ xlc/xlC Compilers
if grep 'xlc' build-aux/depcomp >/dev/null 2>&1
then :
else
  patch -b -p0 2>/dev/null <<\EOF
--- build-aux/depcomp.orig
+++ build-aux/depcomp
@@ -102,6 +102,12 @@
    depmode=msvc7
 fi
 
+if test "$depmode" = xlc; then
+   # IBM C/C++ Compilers xlc/xlC can output gcc-like dependency informations.
+   gccflag=-qmakedep=gcc,-MF
+   depmode=gcc
+fi
+
 case "$depmode" in
 gcc3)
 ## gcc 3 implements dependency tracking that does exactly what
@@ -226,6 +232,13 @@
   rm -f "$tmpdepfile"
   ;;
 
+xlc)
+  # This case exists only to let depend.m4 do its work.  It works by
+  # looking at the text of this script.  This case will never be run,
+  # since it is checked for above.
+  exit 1
+  ;;
+
 aix)
   # The C for AIX Compiler uses -M and outputs the dependencies
   # in a .u file.  In older versions, this file always lives in the
EOF
fi
