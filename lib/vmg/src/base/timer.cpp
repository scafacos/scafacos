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
 * @file   timer.cpp
 * @author Julian Iseringhausen <isering@ins.uni-bonn.de>
 * @date   Tue Sep 6 16:17:40 2011
 *
 * @brief  Class to measure timings.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#ifdef HAVE_MARMOT
#include <enhancempicalls.h>
#include <sourceinfompicalls.h>
#endif
#endif

#include <iostream>
#include <limits>
#include <sstream>

#include "base/helper.hpp"
#include "base/timer.hpp"
#include "comm/comm.hpp"
#include "thirdparty/pugixml/pugixml.hpp"
#include "mg.hpp"

using namespace VMG;

std::map<std::string, TimerData> Timer::td;

void Timer::Start(std::string event)
{
#ifdef HAVE_MPI
#ifdef DEBUG_MEASURE_TIME
  std::map<std::string, TimerData>::iterator iter = td.find(event);
  if (iter == td.end())
    iter = td.insert(std::make_pair(event, TimerData())).first;

  iter->second.time_start = MPI_Wtime();
#endif
#endif
}

void Timer::Stop(std::string event)
{
#ifdef HAVE_MPI
#ifdef DEBUG_MEASURE_TIME
  double time_end = MPI_Wtime();

  std::map<std::string, TimerData>::iterator iter = td.find(event);

  if (time_end - iter->second.time_start < std::numeric_limits<double>::min())
    ++(iter->second.warning);

  ++(iter->second.total);

  iter->second.duration += time_end - iter->second.time_start;
#endif
#endif
}

void Timer::Clear()
{
  td.clear();
}

pugi::xml_node Timer::ToXMLNode()
{
 std::map<std::string, TimerData>::iterator iter;

 pugi::xml_node node_process;
 node_process.append_attribute("Rank").set_value(MG::GetComm()->GlobalRank());

 pugi::xml_node node_timings = node_process.append_child("Timings");

 for (iter=Timer::td.begin(); iter!=Timer::td.end(); ++iter) {

   pugi::xml_node node_entry = node_timings.append_child("Sample");
   node_entry.append_attribute("Name").set_value(Helper::ToString(iter->first).c_str());

   node_entry.append_child("Duration")
     .append_child(pugi::node_pcdata)
     .set_value(Helper::ToString(iter->second.duration).c_str());

   node_entry.append_child("Warnings")
     .append_child(pugi::node_pcdata)
     .set_value(Helper::ToString(iter->second.warning).c_str());

   node_entry.append_child("Total")
     .append_child(pugi::node_pcdata)
     .set_value(Helper::ToString(iter->second.total).c_str());

 }

 return node_process;
}

std::string Timer::ToString()
{
  pugi::xml_node node = Timer::ToXMLNode();
  std::stringstream str;
  node.print(str);
  return str.str();
}

void Timer::Print()
{
#ifdef DEBUG_MEASURE_TIME
 std::map<std::string, TimerData>::const_iterator iter;
 Comm& comm = *MG::GetComm();

 if (comm.GlobalRank() == 0) {
   comm.PrintStringOnce("Running times:");
   for (iter=Timer::td.begin(); iter!=Timer::td.end(); ++iter)
     comm.PrintStringOnce("  %s: %e s (%d)", iter->first.c_str(), iter->second.duration, iter->second.total);
 }
#endif
}

template <class T>
static T min(T* data, int num_data, int& at_rank)
{
  at_rank = 0;
  T min = data[0];

  for (int i=1; i<num_data; ++i)
    if (data[i] < min) {
      at_rank = i;
      min = data[i];
    }
  return min;
}

template <class T>
static T max(T* data, int num_data, int& at_rank)
{
  at_rank = 0;
  T max = data[0];

  for (int i=1; i<num_data; ++i)
    if (data[i] > max) {
      at_rank = i;
      max = data[i];
    }
  return max;
}

template <class T>
static vmg_float avg(T* data, int num_data)
{
  vmg_float average = 0.0;
  vmg_float num_data_inv = 1.0 / static_cast<vmg_float>(num_data);
  for (int i=0; i<num_data; ++i)
    average += data[i] * num_data_inv;
  return average;
}

void Timer::PrintGlobal()
{
#ifdef DEBUG_MEASURE_TIME
  std::map<std::string, TimerData>::const_iterator iter;
  Comm& comm = *MG::GetComm();
  char name[80];

  int rank = comm.GlobalRank();
  int size = comm.GlobalSize();

  vmg_float times[size];
  int calls[size];

  comm.PrintStringOnce("Running times (global):");

  int timer_size = Timer::td.size();
  comm.GlobalBroadcast(timer_size);

  if (rank == 0) {
      for (iter=Timer::td.begin(); iter!=Timer::td.end(); ++iter) {
	std::strcpy(name, iter->first.c_str());
	comm.GlobalBroadcast(name);
	comm.GlobalGather(Timer::td[name].duration, times);
	comm.GlobalGather(Timer::td[name].total, calls);

	int min_calls, max_calls;
	vmg_float avg_calls;
	vmg_float min_duration, max_duration, avg_duration;
	int rank_min_calls, rank_max_calls, rank_min_duration, rank_max_duration;

	min_duration = min(times, size, rank_min_duration);
	max_duration = max(times, size, rank_max_duration);
	avg_duration = avg(times, size);
	min_calls = min(calls, size, rank_min_calls);
	max_calls = max(calls, size, rank_max_calls);
	avg_calls = avg(calls, size);

	comm.PrintStringOnce("  %s: %e s (%d)", iter->first.c_str(), iter->second.duration, iter->second.total);
	comm.PrintStringOnce("    Min: %e s @ %d", min_duration, rank_min_duration);
	comm.PrintStringOnce("    Max: %e s @ %d", max_duration, rank_max_duration);
	comm.PrintStringOnce("    Avg: %e s", avg_duration);
	comm.PrintStringOnce("    Min calls: %d @ %d", min_calls, rank_min_calls);
	comm.PrintStringOnce("    Max calls: %d @ %d", max_calls, rank_max_calls);
	comm.PrintStringOnce("    Avg calls: %f", avg_calls);
    }
  }else {
    for (int i=0; i<timer_size; ++i) {
      comm.GlobalBroadcast(name);
      comm.GlobalGather(Timer::td[name].duration, times);
      comm.GlobalGather(Timer::td[name].total, calls);
    }
  }

#endif
}

std::ostream& VMG::operator<<(std::ostream& out, const Timer&)
{
  return out << Timer::ToString();
}
