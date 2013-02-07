#!/bin/bash

## FFTW has some more macros in m4 directory, while scafacos uses a symbolic link
## Idea: Do not use the symbolic link for FFTWs m4 directory?

FFTW_VERS=3.3
FFTW_DIR=fftw-$FFTW_VERS


# get all files name that are not in .svn directories 
#filelist=$(find $FFTW_DIR -type f | grep -v ".svn")

## change name of libfftw3.a to libfcs_fftw3.a
filelist=$(grep -rl "libfftw3" $FFTW_DIR | grep -v ".svn")
echo "###############################################"
echo "# Replace 'libfftw' prefix in following files #"
echo "###############################################"
echo "filelist =" $filelist
for filename in $filelist; do
  sed -i 's/libfftw3/libfcs_fftw3/g' $filename
done

## delete FFTWs MAINTAINER_MODE build rules, since they use gnu pattern rules %
filelist=$(grep -rl "%\.c" $FFTW_DIR | grep -v ".svn")
echo "#########################################################"
echo "# Delete MAINTAINER_MODE build rules in following files #"
echo "#########################################################"
echo "filelist =" $filelist
for filename in $filelist; do
  sed -i '/if MAINTAINER_MODE/,/endif # MAINTAINER_MODE/d' $filename
done


## do not create fortran wrappers, since we do not need them in ScaFaCoS
sed -i 's/MAINTAINER_MODE/FALSE/g' $FFTW_DIR/mpi/Makefile.am




## api/fftw3.h: change prefix from 'fftw' into 'fsc_fftw'
## no need to change mpi/fftw3-mpi.h
filelist=$(grep -rl "CONCAT(fftw" $FFTW_DIR | grep -v ".svn")
echo "############################################"
echo "# Replace 'fftw' prefix in following files #"
echo "############################################"
echo "filelist =" $filelist
for filename in $filelist; do
  sed -i 's/CONCAT(fftw/CONCAT(fcs_fftw/g' $filename
done

## Enable Maintainer mode per default
## This line does not work correctly
#sed -i "s/AM_MAINTAINER_MODE/AM_MAINTAINER_MODE([enable])/g" $FFTW_DIR/configure.ac

# Add banner to FFTW configure.
sed -i "/AC_INIT/ a\\
\\
AC_MSG_NOTICE([****************************************************************])\\
AC_MSG_NOTICE([*      Configuring FFTW in common/fftw-3.3                     *])\\
AC_MSG_NOTICE([****************************************************************])\\
" $FFTW_DIR/configure.ac

# Add FFTW_MANGLE_PREFIX macro (adds fcs_ to any name).
sed -i "s/\(#define FFTW_MANGLE_QUAD.*\)/\1\n\n#define FFTW_MANGLE_PREFIX(name) FFTW_CONCAT(fcs_, name)/" $FFTW_DIR/api/fftw3.h

# fix inconsistence use of prefixes in FFTW
#sed -i "s/fftw_plan/X(plan)/g" $FFTW_DIR/mpi/f03-wrap.c
#sed -i "s/fftw_r2r_kind/X(r2r_kind)/g" $FFTW_DIR/mpi/f03-wrap.c
#sed -i "s/fftw_complex/X(complex)/g" $FFTW_DIR/mpi/f03-wrap.c
sed -i "s/fftw_plan/fcs_fftw_plan/g" $FFTW_DIR/mpi/f03-wrap.c
sed -i "s/fftw_r2r_kind/fcs_fftw_r2r_kind/g" $FFTW_DIR/mpi/f03-wrap.c
sed -i "s/fftw_complex/fcs_fftw_complex/g" $FFTW_DIR/mpi/f03-wrap.c

sed -i "s/fftw_complex/fcs_fftw_complex/" $FFTW_DIR/api/fftw3.h
sed -i "s/fftwf_complex/fcs_fftwf_complex/" $FFTW_DIR/api/fftw3.h
sed -i "s/fftwl_complex/fcs_fftwl_complex/" $FFTW_DIR/api/fftw3.h
sed -i "s/fftwq_complex/fcs_fftwq_complex/" $FFTW_DIR/api/fftw3.h

sed -i "s/fftw_complex/fcs_fftw_complex/" $FFTW_DIR/mpi/fftw3-mpi.h
sed -i "s/fftwf_complex/fcs_fftwf_complex/" $FFTW_DIR/mpi/fftw3-mpi.h
sed -i "s/fftwl_complex/fcs_fftwl_complex/" $FFTW_DIR/mpi/fftw3-mpi.h

sed -i "s/\(fftw[fl]\?\)_mkstride/fcs_\1_mkstride/" $FFTW_DIR/kernel/ifftw.h
sed -i "s/\(fftw[fl]\?\)_stride_destroy/fcs_\1_stride_destroy/" $FFTW_DIR/kernel/ifftw.h


echo "########################################################"
echo "# Comment installation of pkgconfig in following file: #"
echo "########################################################"
echo "file =" $FFTW_DIR"/Makefile.am"
sed -i "s/pkgconfig/# pkgconfig/" $FFTW_DIR/Makefile.am


## rename check-local target to avoid FFTW test during make check
filelist=$(grep -rl "check-local" $FFTW_DIR | grep -v ".svn")
echo "############################################"
echo "# Replace 'check-local' in following files #"
echo "############################################"
echo "filelist =" $filelist
for filename in $filelist; do
  sed -i "s/check-local/normal-check/g" $filename
done



echo "###########################################"
echo "# Set symbolic links to m4 and build-aux: #"
echo "###########################################"
cd $FFTW_DIR
ln -s ../../../m4 m4
ln -s ../../../build-aux build-aux
cd -

sed -i "/AM_INIT_AUTOMAKE/ iAC_CONFIG_AUX_DIR([build-aux])\\
" $FFTW_DIR/configure.ac

sed -i "s/threads libbench2 . tests mpi doc tools m4/threads libbench2 . tests mpi doc tools/" $FFTW_DIR/Makefile.am
sed -i 's\m4/Makefile\\' $FFTW_DIR/configure.ac


echo "####################################"
echo "# Disable build of docs and tools: #"
echo "####################################"
sed -i -e "/doc\//,/tools\// d" -e "/tools\// d" $FFTW_DIR/configure.ac
sed -i "s/ doc tools//" $FFTW_DIR/Makefile.am


echo "###############################"
echo "# Disable install of headers: #"
echo "###############################"
sed -i -e "s/^\([^#]*include_HEADERS.*\)$/#\1/" -e "s/^\(BUILT_SOURCES.*\)$/#\1/" $FFTW_DIR/api/Makefile.am $FFTW_DIR/mpi/Makefile.am


echo "#########################"
echo "# fix 'make distclean': #"
echo "#########################"
sed -i "/^BUILT_SOURCES/ aDISTCLEANFILES= \$(CODLIST)" \
  $FFTW_DIR/rdft/scalar/r2cb/Makefile.am $FFTW_DIR/rdft/scalar/r2cf/Makefile.am $FFTW_DIR/rdft/scalar/r2r/Makefile.am \
  $FFTW_DIR/rdft/simd/common/Makefile.am \
  $FFTW_DIR/dft/simd/common/Makefile.am \
  $FFTW_DIR/dft/scalar/codelets/Makefile.am

sed -i "/^endif/ a\\
\\
DISTCLEANFILES= \$(EXTRA_DIST)" \
  $FFTW_DIR/rdft/simd/altivec/Makefile.am $FFTW_DIR/rdft/simd/sse2/Makefile.am $FFTW_DIR/rdft/simd/avx/Makefile.am $FFTW_DIR/rdft/simd/neon/Makefile.am \
  $FFTW_DIR/dft/simd/altivec/Makefile.am $FFTW_DIR/dft/simd/sse2/Makefile.am $FFTW_DIR/dft/simd/avx/Makefile.am $FFTW_DIR/dft/simd/neon/Makefile.am


## set svn:ignore to ignore tmp configury files
## newline is necessary after each item
svn propset -R svn:ignore "
Makefile.in
Makefile
" $FFTW_DIR

svn propset svn:ignore "
*.o
*.a
*.mod
.deps
aclocal.m4
autom4te.cache
configure
config.status
config.cache
config.log
config.h.in
config.h.in~
config.h
stamp-*
Makefile.in
Makefile
build-aux
fconfig.h
" $FFTW_DIR

# svn propset svn:ignore "
# Makefile.in
# Makefile
# libtool.m4
# lt~obsolete.m4
# ltsugar.m4
# ltversion.m4
# ltoptions.m4
# " $FFTW_DIR/m4

