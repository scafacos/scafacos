/*
 *    vmg - a versatile multigrid solver
 *    Copyright (C) 2012 Institute for Numerical Simulation, University of Bonn
 *
 *  vmg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  vmg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file   has_request_vec.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Mon Nov 21 13:27:22 2011
 *
 * @brief  A convenience base class for classes using MPI_Requests.
 *
 */


#ifndef HAS_REQUEST_VEC_HPP_
#define HAS_REQUEST_VEC_HPP_

#include <vector>

#ifdef HAVE_MPI

namespace VMG
{

class HasRequestVec
{
protected:

MPI_Request& Request()
{
  request_vec.push_back(MPI_Request());
  return request_vec.back();
}

void WaitAll()
{
  if (!request_vec.empty()) {
#ifndef NDEBUG
    int rval = MPI_Waitall(static_cast<int>(request_vec.size()), &request_vec.front(), MPI_STATUSES_IGNORE);
    assert(rval == MPI_SUCCESS);
#else
    MPI_Waitall(static_cast<int>(request_vec.size()), &request_vec.front(), MPI_STATUSES_IGNORE);
#endif
    request_vec.clear();
  }
}

int TestAll()
{
  int flag = 1;
  if (!request_vec.empty())
    MPI_Testall(request_vec.size(), &request_vec.front(), &flag, MPI_STATUSES_IGNORE);
  return flag;
}

size_t RequestsPending()
{
  return request_vec.size();
}

private:
  std::vector<MPI_Request> request_vec;
};

}

#endif /* HAVE_MPI */

#endif /* HAS_REQUEST_VEC_HPP_ */
