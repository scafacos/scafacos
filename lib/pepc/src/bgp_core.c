/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
* 
* Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
*                         Forschungszentrum Juelich GmbH,
*                         Germany
* 
* PEPC is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* PEPC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
*/

//////////////// BGP-Core Identification //////////////////////
#if defined(__TOS_BGP__)
// see "Using the IBM XL Compilers for Blue Gene", pg. 7:
// "
// __TOS_BGP__ Indicates that the target architecture is PowerPC 450
// "

#include <spi/kernel_interface.h>

int get_my_core()
{
  return Kernel_PhysicalProcessorID();
}

int set_prefetching()
{
  return 0;
}


#elif defined(__TOS_BGQ__)
// see XL C/C++ for Blue Gene/Q, V12.1 > Compiler Reference > 
//       Compiler predefined macros > Macros related to the platform

#include <spi/include/kernel/location.h>
#include <spi/include/l1p/sprefetch.h>

int get_my_core()
{
  return Kernel_ProcessorCoreID();
}

// see http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUQUEEN/UserInfo/CompilingTuning.html#doc1168282bodyText6
int set_prefetching()
{
  return L1P_SetStreamPolicy(L1P_stream_confirmed);
}

#else /* !defined(__TOS_BGP__) && !defined(__TOS_BGQ__) */
int get_my_core()
{
  // we have to be sure that on machines with a standard scheduler
  // no thread thinks, he is sharing its processor with someone else
  static int lastreq = 144;

  return ++lastreq;
}

int set_prefetching()
{
  return 0;
}

#endif

