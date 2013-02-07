/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
* 
* Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
#ifdef __TOS_BGP__
// see "Using the IBM XL Compilers for Blue Gene", pg. 7:
// "
// __TOS_BGP__ Indicates that the target architecture is PowerPC 450
// "

#include <spi/kernel_interface.h>

int get_my_core()
{
  return Kernel_PhysicalProcessorID();
}
#else
int get_my_core()
{
  // we have to be sure that on machines with a standard scheduler
  // no thread thinks, he is sharing its processor with someone else
  static int lastreq = 144;

  return ++lastreq;
}
#endif

