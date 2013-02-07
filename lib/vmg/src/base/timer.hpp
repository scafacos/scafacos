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
 * @file   timer.hpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Sep 6 16:17:40 2011
 *
 * @brief  Class to measure timings.
 *
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <map>
#include <string>

#include "thirdparty/pugixml/pugixml.hpp"

namespace VMG
{

struct TimerData
{
  TimerData():
    time_start(0.0),
    duration(0.0),
    warning(0),
    total(0)
  {}

  double time_start;
  double duration;
  int warning;
  int total;
};

class Timer
{
public:
  Timer() {}
  ~Timer() {}

  static void Start(std::string event);
  static void Stop(std::string event);
  static void Clear();

  static pugi::xml_node ToXMLNode();
  static std::string ToString();
  static void Print();
  static void PrintGlobal();

  static std::map<std::string, TimerData> td;

private:
  static std::string active_timer;
};

std::ostream& operator<<(std::ostream& out, const Timer&);

}

#endif /* TIMER_HPP_ */

