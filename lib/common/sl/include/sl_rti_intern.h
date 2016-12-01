/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Michael Hofmann
 *  
 *  This file is part of ScaFaCoS.
 *  
 *  ScaFaCoS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ScaFaCoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  

 *  
 *  SL - Sorting Library, michael <dot> hofmann <at> informatik <dot> tu-chemnitz <dot> de
 */


#ifndef __SL_RTI_INTERN_H__
#define __SL_RTI_INTERN_H__


/* sl_macro SL_USE_RTI SL_USE_RTI_CMC SL_USE_RTI_TIM SL_USE_RTI_MEM */

#ifndef SL_USE_RTI

 #undef SL_USE_RTI_CMC  /* compare-move-counter */
 #undef SL_USE_RTI_TIM  /* timing */
 #undef SL_USE_RTI_MEM  /* memory */

#endif

#ifdef SL_USE_RTI_CMC

 /* regular commands */
 #define rti_cadd_cmp(n)          (SL_DEFCON(rti).cmc.cmp += n)  /* sl_macro */
 #define rti_cadd_movek(n)        (SL_DEFCON(rti).cmc.movek += n)  /* sl_macro */
 #define rti_cadd_moved(n)        (SL_DEFCON(rti).cmc.moved += n)  /* sl_macro */
 #define rti_cclear_cmp()         (SL_DEFCON(rti).cmc.cmp = 0)  /* sl_macro */
 #define rti_cclear_movek()       (SL_DEFCON(rti).cmc.movek = 0)  /* sl_macro */
 #define rti_cclear_moved()       (SL_DEFCON(rti).cmc.moved = 0)  /* sl_macro */
 #define rti_cclear_all()         (SL_DEFCON(rti).cmc.cmp = SL_DEFCON(rti).cmc.movek = SL_DEFCON(rti).cmc.moved = 0)  /* sl_macro */
 #define rti_ccmp()               my_rti_ccmp(SL_DEFCON(rti))  /* sl_macro */
 #define rti_cmovek()             my_rti_cmovek(SL_DEFCON(rti))  /* sl_macro */
 #define rti_cmoved()             my_rti_cmoved(SL_DEFCON(rti))  /* sl_macro */

 /* chained commands */
 #define cc_rti_cadd_cmp(n)       rti_cadd_cmp(n),  /* sl_macro */
 #define cc_rti_cadd_movek(n)     rti_cadd_movek(n),  /* sl_macro */
 #define cc_rti_cadd_moved(n)     rti_cadd_moved(n),  /* sl_macro */

#else /* SL_USE_RTI_CMC */

 /* regular commands */
 #define rti_cadd_cmp(n)
 #define rti_cadd_movek(n)
 #define rti_cadd_moved(n)
 #define rti_cclear_cmp()
 #define rti_cclear_movek()
 #define rti_cclear_moved()
 #define rti_cclear_all()
 #define rti_ccmp()               0
 #define rti_cmovek()             0
 #define rti_cmoved()             0

 /* chained commands */
 #define cc_rti_cadd_cmp(n)
 #define cc_rti_cadd_movek(n)
 #define cc_rti_cadd_moved(n)

#endif /* SL_USE_RTI_CMC */


#ifdef SL_USE_RTI_TIM

 #define rti_tstart(t)            (SL_DEFCON(rti).tim[t].start = z_time_get_s(), ++SL_DEFCON(rti).tim[t].num)  /* sl_macro */
 #define rti_tstop(t)             (SL_DEFCON(rti).tim[t].stop = z_time_get_s(), SL_DEFCON(rti).tim[t].cumu += (SL_DEFCON(rti).tim[t].last = SL_DEFCON(rti).tim[t].stop - SL_DEFCON(rti).tim[t].start))  /* sl_macro */
 #define rti_tclear(t)            (SL_DEFCON(rti).tim[t].last = 0)  /* sl_macro */
 #define rti_treset(t)            (SL_DEFCON(rti).tim[t].last = SL_DEFCON(rti).tim[t].cumu = 0, SL_DEFCON(rti).tim[t].num = 0)  /* sl_macro */
 #define rti_tlast(t)             my_rti_tlast(SL_DEFCON(rti), t)  /* sl_macro */
 #define rti_tcumu(t)             my_rti_tcumu(SL_DEFCON(rti), t)  /* sl_macro */
 #define rti_tnum(t)              my_rti_tnum(SL_DEFCON(rti), t)  /* sl_macro */

#else

 #define rti_tstart(t)            Z_NOP()
 #define rti_tstop(t)             Z_NOP()
 #define rti_tclear(t)            Z_NOP()
 #define rti_treset(t)            Z_NOP()
 #define rti_tlast(t)             0
 #define rti_tcumu(t)             0
 #define rti_tnum(t)              0

#endif


#ifdef SL_USE_RTI_MEM

 #define rti_minc_alloc()         ++SL_DEFCON(rti).mem.nalloc  /* sl_macro */
 #define rti_minc_free()          ++SL_DEFCON(rti).mem.nfree  /* sl_macro */
 #define rti_malloc(_s_)          (SL_DEFCON(rti).mem.max = z_max(_s_, SL_DEFCON(rti).mem.max), SL_DEFCON(rti).mem.cur += _s_, SL_DEFCON(rti).mem.cur_max = z_max(SL_DEFCON(rti).mem.cur, SL_DEFCON(rti).mem.cur_max))  /* sl_macro */
 #define rti_mfree(_s_)           (SL_DEFCON(rti).mem.cur -= _s_)  /* sl_macro */

 #define cc_rti_minc_alloc()      rti_minc_alloc(),  /* sl_macro */
 #define cc_rti_minc_free()       rti_minc_free(),  /* sl_macro */
 #define cc_rti_malloc(_s_)       rti_malloc(_s_),  /* sl_macro */
 #define cc_rti_mfree(_s_)        rti_mfree(_s_),  /* sl_macro */

#else

 #define rti_minc_alloc()         Z_NOP()
 #define rti_minc_free()          Z_NOP()
 #define rti_malloc(_s_)          Z_NOP()
 #define rti_mfree(_s_)           Z_NOP()

 #define cc_rti_minc_alloc()
 #define cc_rti_minc_free()
 #define cc_rti_malloc(_s_)
 #define cc_rti_mfree(_s_)

#endif


#ifdef SL_USE_RTI
 #define rti_reset()              my_rti_reset(SL_DEFCON(rti))  /* sl_macro */
#else
 #define rti_reset()              Z_NOP()
#endif


#endif /* __SL_RTI_INTERN_H__ */
