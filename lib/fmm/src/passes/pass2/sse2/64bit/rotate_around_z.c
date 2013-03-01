/*
 * Copyright (C) 2012 Ivo Kabadshow, Holger Dachsel
 *
 * This file is part of ScaFaCoS.
 *
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *     
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <xmmintrin.h>

void rotate_around_z(long long nmultipoles,
                     double *romega11, double *romega12, double *romega13, double *romega14,
                     double *iomega11, double *iomega12, double *iomega13, double *iomega14,
                     double *romega21, double *romega22, double *romega23, double *romega24,
                     double *iomega21, double *iomega22, double *iomega23, double *iomega24,
                     double *csmphi, double *csmphipi,
                     double *scr1, double *scr2,
                     double *d2)
{
  int i,j,k,l,m;
  
  __m128d reg00,reg01,reg02,reg03;
  __m128d reg04,reg05,reg06,reg07;
  __m128d reg08,reg09,reg10,reg11;
  __m128d reg12,reg13,reg14,reg15;
  __m128d reg16,reg17;                   /* register for csmphi */
  __m128d reg18,reg19;                   /* register for csmphipi */

  scr1[ 0] = romega11[0];
  scr1[ 1] = romega12[0];
  scr1[ 2] = romega13[0];
  scr1[ 3] = romega14[0];

  scr1[ 4] = iomega11[0];
  scr1[ 5] = iomega12[0];
  scr1[ 6] = iomega13[0];
  scr1[ 7] = iomega14[0];

  scr1[ 8] = romega21[0];
  scr1[ 9] = romega22[0];
  scr1[10] = romega23[0];
  scr1[11] = romega24[0];

  scr1[12] = iomega21[0];
  scr1[13] = iomega22[0];
  scr1[14] = iomega23[0];
  scr1[15] = iomega24[0];

  i = 0;
  j = 0;

  for(l=1;l<=nmultipoles;++l)
  {
    i += 1;
    j += 16;
    k = 2*l;

    _mm_store_pd(&scr1[j   ],_mm_set_pd(romega12[i],romega11[i]));
    _mm_store_pd(&scr1[j+ 2],_mm_set_pd(romega14[i],romega13[i]));
    
    _mm_store_pd(&scr1[j+ 4],_mm_set_pd(iomega12[i],iomega11[i]));
    _mm_store_pd(&scr1[j+ 6],_mm_set_pd(iomega14[i],iomega13[i]));

    _mm_store_pd(&scr1[j+ 8],_mm_set_pd(romega22[i],romega21[i]));
    _mm_store_pd(&scr1[j+10],_mm_set_pd(romega24[i],romega23[i]));

    _mm_store_pd(&scr1[j+12],_mm_set_pd(iomega22[i],iomega21[i]));
    _mm_store_pd(&scr1[j+14],_mm_set_pd(iomega24[i],iomega23[i]));

/*    
    // TODO use mm_set_pd to load/store scr1 
    scr1[j- 1] = romega11[i-1];
    scr1[j   ] = romega12[i-1];
    scr1[j+ 1] = romega13[i-1];
    scr1[j+ 2] = romega14[i-1];

    scr1[j+ 3] = iomega11[i-1];
    scr1[j+ 4] = iomega12[i-1];
    scr1[j+ 5] = iomega13[i-1];
    scr1[j+ 6] = iomega14[i-1];

    scr1[j+ 7] = romega21[i-1];
    scr1[j+ 8] = romega22[i-1];
    scr1[j+ 9] = romega23[i-1];
    scr1[j+10] = romega24[i-1];

    scr1[j+11] = iomega21[i-1];
    scr1[j+12] = iomega22[i-1];
    scr1[j+13] = iomega23[i-1];
    scr1[j+14] = iomega24[i-1];
*/
    for(m=1;m<=k;m+=2)
    {
      i += 1;
      j += 16;

      reg16 = _mm_load_pd(&csmphi[m-1]);     /* load csmphi[m] and csmphi[m+1] into reg16 */

      /* generate [m+1] part in reg17=[csmphi[m+1],csmphi[m+1]]
       * generate reg17=[a1,a1] from [a0,a1] and [a0,a1]  with mask1 = 3*/
      reg17 = _mm_shuffle_pd(reg16, reg16, 0x3);

      /* generate [m] part in reg17=[csmphi[m],csmphi[m]]
       * generate reg16=[a0,a0] from [a0,a1] and [a0,a1]  with mask2*/
      reg16 = _mm_shuffle_pd(reg16, reg16, 0x0); 

      /* load csmphipi[m] and csmphipi[m+1] into reg18 */
      reg18 = _mm_load_pd(&csmphipi[m-1]);   

      /* generate [m+1] part in reg19=[csmphipi[m+1],csmphipi[m+1]]
       * generate reg19=[a1,a1] from [a0,a1] and [a0,a1]  with mask1 = 3*/
      reg19 = _mm_shuffle_pd(reg18, reg18, 0x3);

      /* generate [m] part in reg18=[csmphipi[m],csmphipi[m]]
       * generate reg18=[a0,a0] from [a0,a1] and [a0,a1]  with mask2*/
      reg18= _mm_shuffle_pd(reg18, reg18, 0x0);

      reg00 = _mm_set_pd(romega12[i],romega11[i]);
      reg01 = _mm_set_pd(romega14[i],romega13[i]);

      reg02 = _mm_set_pd(iomega12[i],iomega11[i]);
      reg03 = _mm_set_pd(iomega14[i],iomega13[i]);

      reg04 = _mm_set_pd(romega22[i],romega21[i]);
      reg05 = _mm_set_pd(romega24[i],romega23[i]);

      reg06 = _mm_set_pd(iomega22[i],iomega21[i]);
      reg07 = _mm_set_pd(iomega24[i],iomega23[i]);

      reg08 = _mm_mul_pd(reg00,reg16);
      reg08 = _mm_sub_pd(reg08,_mm_mul_pd(reg02,reg17));

      reg09 = _mm_mul_pd(reg01,reg16);
      reg09 = _mm_sub_pd(reg09,_mm_mul_pd(reg03,reg17));

      reg10 = _mm_mul_pd(reg02,reg16);
      reg10 = _mm_add_pd(reg10,_mm_mul_pd(reg00,reg17));

      reg11 = _mm_mul_pd(reg03,reg16);
      reg11 = _mm_add_pd(reg11,_mm_mul_pd(reg01,reg17));

      reg12 = _mm_mul_pd(reg04,reg18);
      reg12 = _mm_sub_pd(reg12,_mm_mul_pd(reg06,reg19));

      reg13 = _mm_mul_pd(reg05,reg18);
      reg13 = _mm_sub_pd(reg13,_mm_mul_pd(reg07,reg19));

      reg14 = _mm_mul_pd(reg06,reg18);
      reg14 = _mm_add_pd(reg14,_mm_mul_pd(reg04,reg19));

      reg15 = _mm_mul_pd(reg07,reg18);
      reg15 = _mm_add_pd(reg15,_mm_mul_pd(reg05,reg19));

      _mm_store_pd(&scr1[j   ],reg08);
      _mm_store_pd(&scr1[j+ 2],reg09);
      _mm_store_pd(&scr1[j+ 4],reg10);
      _mm_store_pd(&scr1[j+ 6],reg11);
      _mm_store_pd(&scr1[j+ 8],reg12);
      _mm_store_pd(&scr1[j+10],reg13);
      _mm_store_pd(&scr1[j+12],reg14);
      _mm_store_pd(&scr1[j+14],reg15);
    }
  }
}
