#!/bin/bash
csplit -f fmmgradtmain -b %01d.f ../fmmgradt.f '%#include%' '/end subroutine fmmgradt/+1' {0}
mv -v fmmgradtmain0.f fmmgradtmain.f

mv -v fmmgradtmain1.f tmp
csplit -f fmmgradtp -b %02d.f tmp '%#if FMM_MAXNMULTIPOLES %' '/#endif/+1' {*}

rm -v tmp
rm -v fmmgradtp51.f

# modify FMM header file fmmgradt main file
sed -i '1s/#include "fmm.h"/#include "..\/..\/fmm.h"/' fmmgradtmain.f

# include FMM header file into splitted sources
for i in {0{0..0},0{1..9},{10..50}}
do
   echo adding include statement into first line of file fmmgradtp${i}.f
   sed -i 1i'#include "../../fmm.h"' fmmgradtp${i}.f
done

