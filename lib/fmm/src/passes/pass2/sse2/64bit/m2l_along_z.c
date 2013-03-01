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

void m2l_along_z(long long nmultipoles, double *scr1, double *scr2, double *d2, double *fr, double *sg)
{
  int mmmm,mmm,mm,m;
  int i,j,k,l,n,nn;

  __m128d reg00,reg01,reg02,reg03;
  __m128d reg04,reg05,reg06,reg07;
  __m128d reg08,reg09,reg10,reg11;
  __m128d reg12,reg13,reg14,reg15;
  __m128d reg16,reg17;                   /* register for rotation matrix TODO: rename regdmat1,regdmat2*/
  __m128d reg18,reg19;                   /* register for g,gl,glm */

  i = -15;

  __m128d regzero = _mm_setzero_pd();

  reg08 = regzero;
  reg09 = regzero;
  reg10 = regzero;
  reg11 = regzero;

  for(j=0;j<=nmultipoles;++j)
  {
    i += 16;

    reg00 = _mm_load_pd(&scr2[i-1]);
    reg01 = _mm_load_pd(&scr2[i+1]);
    reg04 = _mm_load_pd(&scr2[i+7]);
    reg05 = _mm_load_pd(&scr2[i+9]);

    reg18 = _mm_load1_pd(&fr[j]);

    reg08 = _mm_add_pd(reg08,_mm_mul_pd(reg00,reg18));
    reg09 = _mm_add_pd(reg09,_mm_mul_pd(reg01,reg18));
    reg12 = _mm_add_pd(reg12,_mm_mul_pd(reg04,reg18));
    reg13 = _mm_add_pd(reg13,_mm_mul_pd(reg05,reg18));
  }

  _mm_store_pd(&scr1[ 0],reg12);
  _mm_store_pd(&scr1[ 2],reg13);
  _mm_store_pd(&scr1[ 4],regzero);
  _mm_store_pd(&scr1[ 6],regzero);
  _mm_store_pd(&scr1[ 8],reg08);
  _mm_store_pd(&scr1[10],reg09);
  _mm_store_pd(&scr1[12],regzero);
  _mm_store_pd(&scr1[14],regzero);

  i = 1;

  for(l=1;l<=nmultipoles;++l)
  {
    i += 16 * l;
    j = -15;
    k = nmultipoles+l;

    reg08 = regzero;
    reg09 = regzero;
    reg12 = regzero;
    reg13 = regzero;

    for(m=l;m<=k;++m)
    {
      j += 16;

      reg00 = _mm_load_pd(&scr2[j-1]);
      reg01 = _mm_load_pd(&scr2[j+1]);
      reg04 = _mm_load_pd(&scr2[j+7]);
      reg05 = _mm_load_pd(&scr2[j+9]);

      reg18 = _mm_load1_pd(&fr[m]);

      reg08 = _mm_add_pd(reg08,_mm_mul_pd(reg00,reg18));
      reg09 = _mm_add_pd(reg09,_mm_mul_pd(reg01,reg18));
      reg12 = _mm_add_pd(reg12,_mm_mul_pd(reg04,reg18));
      reg13 = _mm_add_pd(reg13,_mm_mul_pd(reg05,reg18));
    }

    reg18 = _mm_load1_pd(&sg[l]);

    reg12 = _mm_mul_pd(reg12,reg18);
    _mm_store_pd(&scr1[i- 1],reg12);

    reg13 = _mm_mul_pd(reg13,reg18);
    _mm_store_pd(&scr1[i+ 1],reg13);

    _mm_store_pd(&scr1[i+ 3],regzero);
    _mm_store_pd(&scr1[i+ 5],regzero);

    reg08 = _mm_mul_pd(reg08,reg18);
    _mm_store_pd(&scr1[i+ 7],reg08);

    reg09 = _mm_mul_pd(reg09,reg18);
    _mm_store_pd(&scr1[i+ 9],reg09);

    _mm_store_pd(&scr1[i+11],regzero);
    _mm_store_pd(&scr1[i+13],regzero);
  }

  mm = 16 * nmultipoles;

  i = 1;
  n = mm+1;

  for(m=1;m<=nmultipoles;++m)
  {
    i += 16 * m;
    j = i;

    for(l=m;l<=nmultipoles;++l)
    {

      j += 16 * l;
      nn = n;
      k = m + l;
      mmm = nmultipoles + l;

      reg08 = regzero;
      reg09 = regzero;
      reg10 = regzero;
      reg11 = regzero;
      reg12 = regzero;
      reg13 = regzero;
      reg14 = regzero;
      reg15 = regzero;

      for(mmmm=k;mmmm<=mmm;++mmmm)
      {
        nn += 16;

        reg00 = _mm_load_pd(&scr2[nn- 1]);
        reg01 = _mm_load_pd(&scr2[nn+ 1]);
        reg02 = _mm_load_pd(&scr2[nn+ 3]);
        reg03 = _mm_load_pd(&scr2[nn+ 5]);
        reg04 = _mm_load_pd(&scr2[nn+ 7]);
        reg05 = _mm_load_pd(&scr2[nn+ 9]);
        reg06 = _mm_load_pd(&scr2[nn+11]);
        reg07 = _mm_load_pd(&scr2[nn+13]);

        reg18 = _mm_load1_pd(&fr[mmmm]);

        reg08 = _mm_add_pd(reg08,_mm_mul_pd(reg00,reg18));
        reg09 = _mm_add_pd(reg09,_mm_mul_pd(reg01,reg18));

        reg10 = _mm_sub_pd(reg10,_mm_mul_pd(reg02,reg18));
        reg11 = _mm_sub_pd(reg11,_mm_mul_pd(reg03,reg18));

        reg12 = _mm_add_pd(reg12,_mm_mul_pd(reg04,reg18));
        reg13 = _mm_add_pd(reg13,_mm_mul_pd(reg05,reg18));

        reg14 = _mm_sub_pd(reg14,_mm_mul_pd(reg06,reg18));
        reg15 = _mm_sub_pd(reg15,_mm_mul_pd(reg07,reg18));
      }

      reg18 = _mm_load1_pd(&sg[k]);

      reg12 = _mm_mul_pd(reg12,reg18);
      _mm_store_pd(&scr1[j- 1],reg12);

      reg13 = _mm_mul_pd(reg13,reg18);
      _mm_store_pd(&scr1[j+ 1],reg13);

      reg14 = _mm_mul_pd(reg14,reg18);
      _mm_store_pd(&scr1[j+ 3],reg14);

      reg15 = _mm_mul_pd(reg15,reg18);
      _mm_store_pd(&scr1[j+ 5],reg15);

      reg08 = _mm_mul_pd(reg08,reg18);
      _mm_store_pd(&scr1[j+ 7],reg08);

      reg09 = _mm_mul_pd(reg09,reg18);
      _mm_store_pd(&scr1[j+ 9],reg09);

      reg10 = _mm_mul_pd(reg10,reg18);
      _mm_store_pd(&scr1[j+11],reg10);

      reg11 = _mm_mul_pd(reg11,reg18);
      _mm_store_pd(&scr1[j+13],reg11);
    }

    n += mm;
    mm -= 16;
  }
}
