/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/* no ghosts and not periodic */
#undef GRIDSORT_FRONT_TPROC_GHOST
#undef GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_bounds
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_bounds_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_zslices_zonly
#include "gridsort_front_tproc.h"

/* with ghosts and not periodic */
#define GRIDSORT_FRONT_TPROC_GHOST
#undef GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_bounds
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_bounds_tricl
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_zonly
#include "gridsort_front_tproc.h"*/

/* no ghosts and periodic */
#undef GRIDSORT_FRONT_TPROC_GHOST
#define GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_bounds
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_bounds_tricl
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#define GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_periodic_zslices_zonly
#include "gridsort_front_tproc.h"*/

/* with ghosts and periodic */
#define GRIDSORT_FRONT_TPROC_GHOST
#define GRIDSORT_FRONT_TPROC_PERIODIC

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_tricl
#include "gridsort_front_tproc.h"

#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_bounds
#include "gridsort_front_tproc.h"

#define GRIDSORT_FRONT_TPROC_TRICLINIC
#define GRIDSORT_FRONT_TPROC_BOUNDS
#undef GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_bounds_tricl
#include "gridsort_front_tproc.h"

/*#undef GRIDSORT_FRONT_TPROC_TRICLINIC
#undef GRIDSORT_FRONT_TPROC_BOUNDS
#define GRIDSORT_FRONT_TPROC_ZSLICES
#undef GRIDSORT_FRONT_TPROC_ZONLY
#define GRIDSORT_FRONT_TPROC_NAME  gridsort_front_tproc_ghost_periodic_zslices
#include "gridsort_front_tproc.h"*/
