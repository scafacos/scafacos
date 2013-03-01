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

void rotate_around_y(long long nmultipoles, double *scr1, double *scr2, double *d2)
{
  long long mmmmmm,mmmmm,mmmm,mmm,mm,m;
  long long i,j,k,l,n,nn,nnn;

  __m128d reg00,reg01,reg02,reg03;
  __m128d reg04,reg05,reg06,reg07;
  __m128d reg08,reg09,reg10,reg11;
  __m128d reg12,reg13,reg14,reg15;
  __m128d reg16,reg17;                   /* register for rotation matrix  TODO rename regdmat1,regdmat2*/
  __m128d reg18,reg19,reg20;             /* register for gl,g,glm */

  const __m128d regzero = _mm_setzero_pd();
  const __m128d regone = _mm_set_pd(1.0e0l,1.0e0l);

  mmmm = nmultipoles + 1;
  mmm = 0;

  _mm_store_pd(&scr2[ 0],_mm_load_pd(&scr1[ 0]));
  _mm_store_pd(&scr2[ 2],_mm_load_pd(&scr1[ 2]));
  _mm_store_pd(&scr2[ 4],_mm_load_pd(&scr1[ 4]));
  _mm_store_pd(&scr2[ 6],_mm_load_pd(&scr1[ 6]));
  _mm_store_pd(&scr2[ 8],_mm_load_pd(&scr1[ 8]));
  _mm_store_pd(&scr2[10],_mm_load_pd(&scr1[10]));
  _mm_store_pd(&scr2[12],_mm_load_pd(&scr1[12]));
  _mm_store_pd(&scr2[14],_mm_load_pd(&scr1[14]));

  i = 1;
  j = 1;
  k = 1;

  reg18 = regone;

  for(l=1;l<=nmultipoles;++l)
  {
    reg18 = _mm_sub_pd(regzero, reg18);  /* negate reg18, gl*/
    i += 1;
    j += 16;
    k += 16 * l;
    n = k;

    mmm += 1;

    reg00 = _mm_load_pd(&scr1[n-1]);
    reg01 = _mm_load_pd(&scr1[n+1]);

    reg04 = _mm_load_pd(&scr1[n+7]);
    reg05 = _mm_load_pd(&scr1[n+9]);
    
    reg17 = _mm_load1_pd(&d2[mmm-1]);   /* loading roation matrix with element mmm-1 twice into reg17 -> TODO: reg17 could be reg16??*/

    reg08 = _mm_mul_pd(reg00,reg17);
    reg09 = _mm_mul_pd(reg01,reg17);
    
    reg12 = _mm_mul_pd(reg04,reg17);
    reg12 = _mm_mul_pd(reg12,reg18);

    reg13 = _mm_mul_pd(reg05,reg17);
    reg13 = _mm_mul_pd(reg13,reg18);

    mmm += 1;                        /* always allocated, when entering this routine */

    nn = n + 16;
    n += 16 * l;
    reg19 = reg18; /* reg19 = g */

    for(mm=nn;mm<=n;mm+=16)
    {
      reg19 = _mm_sub_pd(regzero, reg19);  /* negate reg19*/

      mmm += 1;

      reg00 = _mm_load_pd(&scr1[mm-1]);
      reg01 = _mm_load_pd(&scr1[mm+1]);

      reg04 = _mm_load_pd(&scr1[mm+7]);
      reg05 = _mm_load_pd(&scr1[mm+9]);

      reg17 = _mm_load1_pd(&d2[mmm-1]);   /* loading rotation matrix with element mmm-1 twice into reg17 */

      reg08 = _mm_add_pd(reg08,_mm_mul_pd(reg00,reg17));
      reg09 = _mm_add_pd(reg09,_mm_mul_pd(reg01,reg17));

      reg04 = _mm_mul_pd(reg04,reg19);
      reg12 = _mm_add_pd(reg12,_mm_mul_pd(reg04,reg17));

      reg05 = _mm_mul_pd(reg05,reg19);
      reg13 = _mm_add_pd(reg13,_mm_mul_pd(reg05,reg17));

      mmm += 1;                      /* always allocated, when entering this routine */
    }
      
    _mm_store_pd(&scr2[j- 1],reg08);
    _mm_store_pd(&scr2[j+ 1],reg09);

    _mm_store_pd(&scr2[j+ 3],regzero);
    _mm_store_pd(&scr2[j+ 5],regzero);

    _mm_store_pd(&scr2[j+ 7],reg12);
    _mm_store_pd(&scr2[j+ 9],reg13);
    
    _mm_store_pd(&scr2[j+11],regzero);
    _mm_store_pd(&scr2[j+13],regzero);

    mmmmm = mmmm;
    mmmmmm = i;

    reg20 = _mm_sub_pd(regzero, reg18);  /* reg20 = glm */
    
    for(m=1;m<=l;++m)
    {
      n = k;
      
      reg00 = _mm_load_pd(&scr1[n- 1]);
      reg01 = _mm_load_pd(&scr1[n+ 1]);
      reg02 = _mm_load_pd(&scr1[n+ 3]);
      reg03 = _mm_load_pd(&scr1[n+ 5]);
      reg04 = _mm_load_pd(&scr1[n+ 7]);
      reg05 = _mm_load_pd(&scr1[n+ 9]);
      reg06 = _mm_load_pd(&scr1[n+11]);
      reg07 = _mm_load_pd(&scr1[n+13]);
      
      reg16 = _mm_load_pd(&d2[mmm]);  /* load mmm+1 and mmm+2 into reg16 */
      
      /* generate imag part in reg17=[imag,imag] */
      reg17 = _mm_shuffle_pd(reg16, reg16, 0x3); /* generate reg17=[a1,a1] from [a0,a1] and [a0,a1]  with mask1=3 */

      /* generate real part in reg17=[real, real] */
      reg16 = _mm_shuffle_pd(reg16, reg16, 0x0); /* generate reg16=[a0,a0] from [a0,a1] and [a0,a1]  with mask2=0 */
      
      reg08 = _mm_mul_pd(reg00,reg16);
      reg09 = _mm_mul_pd(reg01,reg16);
      reg10 = _mm_mul_pd(reg02,reg17);
      reg11 = _mm_mul_pd(reg03,reg17);

      reg04 = _mm_mul_pd(reg04,reg20);
      reg12 = _mm_mul_pd(reg04,reg16);

      reg05 = _mm_mul_pd(reg05,reg20);
      reg13 = _mm_mul_pd(reg05,reg16);

      reg19 = reg20; /* generate g */
      reg20 = _mm_sub_pd(regzero, reg20);  /* negate glm to -glm*/

      reg06 = _mm_mul_pd(reg06,reg20);
      reg14 = _mm_mul_pd(reg06,reg17);

      reg07 = _mm_mul_pd(reg07,reg20);
      reg15 = _mm_mul_pd(reg07,reg17);

      mmm += 2;
      nn = n + 16;
      n += 16 * l;

      for(nnn=nn;nnn<=n;nnn+=16)
        {
          reg19 = _mm_sub_pd(regzero, reg19);  /* negate reg19 = g*/
	  
          reg00 = _mm_load_pd(&scr1[nnn- 1]);
          reg01 = _mm_load_pd(&scr1[nnn+ 1]);
          reg02 = _mm_load_pd(&scr1[nnn+ 3]);
          reg03 = _mm_load_pd(&scr1[nnn+ 5]);
          reg04 = _mm_load_pd(&scr1[nnn+ 7]);
          reg05 = _mm_load_pd(&scr1[nnn+ 9]);
          reg06 = _mm_load_pd(&scr1[nnn+11]);
          reg07 = _mm_load_pd(&scr1[nnn+13]);

          reg16 = _mm_load_pd(&d2[mmm]);  /* load mmm+1 and mmm+2 into reg16 */
      
          /* generate imag part in reg17=[imag,imag] */
          reg17 = _mm_shuffle_pd(reg16, reg16, 0x3); /* generate reg17=[a1,a1] from [a0,a1] and [a0,a1]  with mask1=3 */

          /* generate real part in reg17=[real, real] */
          reg16 = _mm_shuffle_pd(reg16, reg16, 0x0); /* generate reg16=[a0,a0] from [a0,a1] and [a0,a1]  with mask2=0 */

          reg08 = _mm_add_pd(reg08,_mm_mul_pd(reg00,reg16));
          reg09 = _mm_add_pd(reg09,_mm_mul_pd(reg01,reg16));

          reg10 = _mm_add_pd(reg10,_mm_mul_pd(reg02,reg17));
          reg11 = _mm_add_pd(reg11,_mm_mul_pd(reg03,reg17));

          reg04 = _mm_mul_pd(reg04,reg19);
          reg12 = _mm_add_pd(reg12,_mm_mul_pd(reg04,reg16));

          reg05 = _mm_mul_pd(reg05,reg19);
          reg13 = _mm_add_pd(reg13,_mm_mul_pd(reg05,reg16));

          reg06 = _mm_mul_pd(reg06,reg19);
          reg14 = _mm_sub_pd(reg14,_mm_mul_pd(reg06,reg17));
	  
          reg07 = _mm_mul_pd(reg07,reg19);
          reg15 = _mm_sub_pd(reg15,_mm_mul_pd(reg07,reg17));

          mmm += 2;

        }

        mmmmmm += mmmmm;
        mmmmm -= 1;
        mm = mmmmmm - m;

        nn = 16 * mm - 15;

        _mm_store_pd(&scr2[nn- 1],reg08);
        _mm_store_pd(&scr2[nn+ 1],reg09);
        _mm_store_pd(&scr2[nn+ 3],reg10);
        _mm_store_pd(&scr2[nn+ 5],reg11);
        _mm_store_pd(&scr2[nn+ 7],reg12);
        _mm_store_pd(&scr2[nn+ 9],reg13);
        _mm_store_pd(&scr2[nn+11],reg14);
        _mm_store_pd(&scr2[nn+13],reg15);

    }
  }
}
