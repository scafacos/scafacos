#!/bin/bash
csplit -f fmmmopnmain -b %01d.f ../fmmmopn.f '%#include%' '/end subroutine fmmmopn/+1' {0}
mv -v fmmmopnmain0.f fmmmopnmain.f

mv -v fmmmopnmain1.f tmp
csplit -f fmmmopnp -b %02d.f tmp '%#if FMM_MAXNMULTIPOLES %' '/#endif/+1' {*}

rm -v tmp
rm -v fmmmopnp51.f

# modify FMM header file fmmmopn main file
sed -i '1s/#include "fmm.h"/#include "..\/..\/fmm.h"/' fmmmopnmain.f

# include FMM header file into splitted sources
for i in {0{0..0},0{1..9},{10..50}}
do
   echo adding include statement into first line of file fmmmopnp${i}.f
   sed -i 1i'#include "../../fmm.h"' fmmmopnp${i}.f
done


