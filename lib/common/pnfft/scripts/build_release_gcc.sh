#!/bin/sh -e

if [ "x$BUILD_RELEASE_CALLED_FROM_MAKE_RELEASE" != "xyes" ]; then
  echo "!!! Error: This is a helper script that should not be called individually. Use the corresponding make_release script instead! !!!"
  exit
fi

myprefix=$HOME/local
PNFFT_VERSION=1.0.5-alpha
PFFT_VERSION=1.0.6-alpha
FFTW_VERSION=3.3.3
INSTDIR=$myprefix/pnfft-$PNFFT_VERSION
PFFTDIR=$myprefix/pfft-$PFFT_VERSION
FFTWDIR=$myprefix/fftw-$FFTW_VERSION
TMP="tmp-pnfft-$PNFFT_VERSION"

# bash check if directory exists
if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
	read answer
	if [ ${answer} = "y" ]; then
		rm -rf $TMP
	else
		echo "Program aborted."
		exit 1
	fi
fi

mkdir $TMP && cd $TMP
cd ../.. && ./bootstrap.sh && cd -
../../configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR --with-pfft=$PFFTDIR --disable-shared

make -j 4 distcheck

