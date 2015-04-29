/*
 *  Copyright (C) 2011, 2012, 2013, 2014, 2015 Michael Hofmann
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


#ifndef __SL_TYPES_INTERN_H__
#define __SL_TYPES_INTERN_H__


#define lb_bin_count(_lb_, _b_, _j_)    ((_lb_)->bins[(_b_) * (_lb_)->nelements + _j_].s.size)
#ifdef elem_weight
# define lb_bin_weight(_lb_, _b_, _j_)  ((_lb_)->bins[(_b_) * (_lb_)->nelements + _j_].weight)
#else
# define lb_bin_weight(_lb_, _b_, _j_)  0
#endif


#if defined(elem_weight) && defined(sl_weight_intequiv)
# define lb_bin_counts(_lb_, _b_, _j_, _k_)    ((_lb_)->bin_cws + z_get4d((_lb_)->bcws[_b_], (_lb_)->bin_cw_factor, 0, (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
# define lb_bin_weights(_lb_, _b_, _j_, _k_)   ((_lb_)->bin_cws + z_get4d((_lb_)->bcws[_b_], (_lb_)->bin_cw_factor, 1, (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
#else
# define lb_bin_counts(_lb_, _b_, _j_, _k_)    ((_lb_)->bin_cs + z_get3d((_lb_)->bcws[_b_], (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
# ifdef elem_weight
#  define lb_bin_weights(_lb_, _b_, _j_, _k_)  ((_lb_)->bin_ws + z_get3d((_lb_)->bcws[_b_], (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
# else
#  define lb_bin_weights(_lb_, _b_, _j_, _k_)  NULL
# endif
#endif


#if defined(elem_weight) && defined(sl_weight_intequiv)
# define lb_counts(_lb_, _b_, _k_)     ((_lb_)->cws        + z_get3d((_lb_)->bcws[_b_], (_lb_)->cw_factor, 0,               (_lb_)->bm->nbins, (_k_)))
# define lb_weights(_lb_, _b_, _k_)    ((_lb_)->cws        + z_get3d((_lb_)->bcws[_b_], (_lb_)->cw_factor, (_lb_)->w_index, (_lb_)->bm->nbins, (_k_)))
# define lb_prefix_count(_lb_, _b_)    ((_lb_)->prefix_cws + z_get2d((_lb_)->bcws[_b_], (_lb_)->cw_factor, 0))
# define lb_prefix_weight(_lb_, _b_)   ((_lb_)->prefix_cws + z_get2d((_lb_)->bcws[_b_], (_lb_)->cw_factor, (_lb_)->w_index))
#else
# define lb_counts(_lb_, _b_, _k_)     ((_lb_)->cs        + z_get2d((_lb_)->bcws[_b_], (_lb_)->bm->nbins, (_k_)))
# define lb_prefix_count(_lb_, _b_)    ((_lb_)->prefix_cs + z_get1d((_lb_)->bcws[_b_]))
# ifdef elem_weight
#  define lb_weights(_lb_, _b_, _k_)   ((_lb_)->ws        + z_get2d((_lb_)->bcws[_b_], (_lb_)->bm->nbins, (_k_)))
#  define lb_prefix_weight(_lb_, _b_)  ((_lb_)->prefix_ws + z_get1d((_lb_)->bcws[_b_]))
# else
#  define lb_weights(_lb_, _b_, _k_)   NULL
#  define lb_prefix_weight(_lb_, _b_)  NULL
# endif
#endif


#if defined(elem_weight) && defined(sl_weight_intequiv)
# define gb_counts(_gb_, _b_, _k_)     ((_gb_)->cws        + z_get3d((_gb_)->bcws[_b_], (_gb_)->lb.cw_factor, 0,                  (_gb_)->bm->nbins, (_k_)))
# define gb_weights(_gb_, _b_, _k_)    ((_gb_)->cws        + z_get3d((_gb_)->bcws[_b_], (_gb_)->lb.cw_factor, (_gb_)->lb.w_index, (_gb_)->bm->nbins, (_k_)))
# define gb_prefix_count(_gb_, _b_)    ((_gb_)->prefix_cws + z_get2d((_gb_)->bcws[_b_], (_gb_)->lb.cw_factor, 0))
# define gb_prefix_weight(_gb_, _b_)   ((_gb_)->prefix_cws + z_get2d((_gb_)->bcws[_b_], (_gb_)->lb.cw_factor, (_gb_)->lb.w_index))
#else
# define gb_counts(_gb_, _b_, _k_)     ((_gb_)->cs        + z_get2d((_gb_)->bcws[_b_], (_gb_)->bm->nbins, (_k_)))
# define gb_prefix_count(_gb_, _b_)    ((_gb_)->prefix_cs + z_get1d((_gb_)->bcws[_b_]))
# ifdef elem_weight
#  define gb_weights(_gb_, _b_, _k_)   ((_gb_)->ws        + z_get2d((_gb_)->bcws[_b_], (_gb_)->bm->nbins, (_k_)))
#  define gb_prefix_weight(_gb_, _b_)  ((_gb_)->prefix_ws + z_get1d((_gb_)->bcws[_b_]))
# else
#  define gb_weights(_gb_, _b_, _k_)   NULL
#  define gb_prefix_weight(_gb_, _b_)  NULL
# endif
#endif


#endif /* __SL_TYPES_INTERN_H__ */
