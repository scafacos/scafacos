#!/bin/bash
csplit -f fmmoopnmain -b %01d.f ../fmmoopn.f '%#include%' '/end subroutine fmmoopn/+1' {0}
mv -v fmmoopnmain0.f fmmoopnmain.f

mv -v fmmoopnmain1.f tmp
csplit -f fmmoopnp -b %02d.f tmp '%#if FMM_MAXNMULTIPOLES %' '/#endif/+1' {*}

rm -v tmp
rm -v fmmoopnp51.f

# modify FMM header file fmmoopn main file
sed -i '1s/#include "fmm.h"/#include "..\/..\/fmm.h"/' fmmoopnmain.f

# include FMM header file into splitted sources
for i in {0{0..0},0{1..9},{10..50}}
do
   echo adding include statement into first line of file fmmoopnp${i}.f
   sed -i 1i'#include "../../fmm.h"' fmmoopnp${i}.f
done

