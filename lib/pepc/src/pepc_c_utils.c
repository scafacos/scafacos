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

/*************************************************************************
>
>  different utils for directly accessing POSIX functions
>
*************************************************************************/

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/stat.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>

#ifndef PATH_MAX 
  #define MYPATH_MAX 255 
#else 
  #define MYPATH_MAX PATH_MAX 
#endif 


void create_directory_c(char dirname[])
{
  char cwd[MYPATH_MAX];
  char fullpath[MYPATH_MAX];

  getcwd(cwd, MYPATH_MAX);
  sprintf(fullpath, "%s/%s", cwd, dirname);

  if (0 != mkdir(fullpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
  {
    // ignore the error if the directory already existed
    if (EEXIST != errno) printf("Error while creating directory %s: %d\n", fullpath, errno);
  }
}


