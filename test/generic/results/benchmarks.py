import numpy
import xml.dom.minidom
import sys
import traceback
import re
import matplotlib.pylab as plt

method2format = {
    'fmm': { 'linestyle': '--', 'color': 'b', 'marker': 'v', 'linewidth': 2},
    'memd': { 'linestyle': '-', 'color': 'g', 'marker': '>', 'linewidth': 2},
    'p2nfft': { 'linestyle': '--', 'color': 'r', 'marker': '^', 'linewidth': 2},
    'p3m': { 'linestyle': '-', 'color': 'c', 'marker': 'o', 'linewidth': 2},
    'pp3mg': { 'linestyle': '--', 'color': 'k', 'marker': '<', 'linewidth': 2},
    'vmg': { 'linestyle': '-', 'color': 'y', 'marker': 'D', 'linewidth': 2},
    'ewald': { 'linestyle': '-', 'color': 'b', 'marker': 'x', 'linewidth': 2},
    'pepc': { 'linestyle': '-', 'color': 'r', 'marker': 'x', 'linewidth': 2},
    }

def fmt(method):
    if (method in method2format): return method2format[method]
    else: return { 'linestyle': '-', 'color': 'g' }

class Data:
    pass
        
def read(filename):
    """Read the XML benchmark file filename.

    Returns a list consisting out of the following elements:
    * a numpy 3D-array of floats with the following axes:
      - the number of charges in the system
      - the tolerance of the benchmark
      - the number of cores
    """

    # Method, Testcase, #charges -> array((#cores, time))
    file = open(filename, "r")
    doc = xml.dom.minidom.parse(file)
    file.close()

    # get top-level element
    benchmarks_el = doc.documentElement

    timing = {}
    speedup = {}
    efficiency = {}

    timing_counter = 0
    dataset_counter = 0

    print "%4s %20s %9s %5s" % ("#", "Method", "#q", "tol.")
    print "%4s %20s %9s %5s" % ("----", "--------", "--", "----")

    all_methods = set()
    all_charges = set()
    all_tolerances = set()
    all_cores = set()

    # loop over all timings
    for timings_el in benchmarks_el.getElementsByTagName("timings"):
        timing_counter += 1
        # get method, charges and tolerance
        try:
            methodname = timings_el.attributes['method'].value.lower()
            charges = int(timings_el.attributes['charges'].value)
            tolerance = float(timings_el.attributes['tolerance'].value)
        except KeyError as error:
            print "Error reading mandatory attribute %s of <timings>-element #%d:" % (error, timings_counter)
            print timings_el.toxml()[:100] + "..."
            print "Skipping <timings> element #%d." % timing_counter
            continue

        all_methods.add(methodname)
        all_charges.add(charges)
        all_tolerances.add(tolerance)

        if methodname not in timing:
            timing[methodname] = {}
        if charges not in timing[methodname]:
            timing[methodname][charges] = {}
        if tolerance not in timing[methodname][charges]:
            current_dict = {}
            timing[methodname][charges][tolerance] = current_dict         
        else:
            print "   Dataset already exists: "
            print "     method=%s charges=%d tolerance=%e" % (methodname, charges, tolerance)
            print "   Skipping <timings> element #%d." % timing_counter
            continue

        # can read dataset
        print "%4d %20s %9s %5s" % (dataset_counter, methodname, charges, tolerance)

        data_els = timings_el.getElementsByTagName("data")
        if len(data_els) < 1:
            print "   Error reading mandatory element <data> of <timings>-element #%d:" % timing_counter
            print "   " + timings_el.toxml()[:100] + "..."
            print "   Skipping <timings> element #%d." % timing_counter
            continue
        elif len(data_els) > 1:
            print "   <timings>-element #%d had more than one <data> child:" % timing_counter
            print "   " + timings_el.toxml()[:100] + "..."
            print "   Skipping <timings> element #%d." % timing_counter
            continue

        datatext = data_els[0].firstChild.data
        lines = datatext.split('\n')
        for line in lines:
            line = line.strip()
            if line == '': continue
            match = re.match('(\d+)[\s,:]+((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)', line)
            if match:
                cores = int(match.group(1))
                time = float(match.group(2))
                all_cores.add(cores)
                current_dict[cores] = time
            else:
                print "   skipping bad data line %s" % line

        dataset_counter += 1

    print "Finished reading data file."
    print "Got %d datasets." % dataset_counter

    all_charges = list(all_charges)
    all_charges.sort()
    all_tolerances = list(all_tolerances)
    all_tolerances.sort(reverse=True)
    all_cores = list(all_cores)
    all_cores.sort()

    # Replace inner dicts with numpy arrays
    for methodname in timing:
        # create arrays
        timing_array = numpy.empty((len(all_charges), len(all_tolerances), len(all_cores)))
        timing_array[:,:,:] = numpy.NaN
        # fill array
        for charges in timing[methodname]:
            ixcha = all_charges.index(charges)
            for tolerance in timing[methodname][charges]:
                ixtol = all_tolerances.index(tolerance)
                for cores in timing[methodname][charges][tolerance]:
                    try:
                      ixcor = all_cores.index(cores)
                      timing_array[ixcha,ixtol,ixcor] = \
                          timing[methodname][charges][tolerance][cores]
                    except ValueError:
                      print "cores = %d not found in list." % cores

        timing[methodname] = timing_array

        speedup_array = timing_array.copy()
        efficiency_array = timing_array.copy()

        # find minimal number of cores with a defined time
        mincoresix = 0
        while mincoresix < timing_array.shape[2]-1 and \
          numpy.isnan(timing_array[ixcha,ixtol,mincoresix]): 
                mincoresix += 1

        # compute estimated 1-core-time
        est_serial_time = \
            timing_array[ixcha,ixtol,mincoresix] * all_cores[mincoresix]

        # speedup
        speedup_array[ixcha,ixtol,:] = est_serial_time
        speedup_array[ixcha,ixtol,:] /= timing_array[ixcha,ixtol,:]

        # efficiency
        efficiency_array[ixcha,ixtol,:] = \
            speedup_array[ixcha,ixtol,:] / all_cores

        speedup[methodname] = speedup_array
        efficiency[methodname] = efficiency_array

    return { 'methods' : all_methods,
             'charges' : all_charges,
             'cores' : all_cores,
             'tolerances' : all_tolerances,
             'timing' : timing,
             'speedup' : speedup,
             'efficiency' : efficiency }
