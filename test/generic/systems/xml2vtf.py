#!/usr/bin/python
#
# Script that transforms the testcase format into the VTF format.
# Note that at the moment it cannot handle cases where the number of
# atoms changes between configurations. Also, the box vectors
# currently have to be parallel to the main axes.
#
import gzip
import sys
import xml.dom.minidom

for filename in sys.argv[1:]:
    # open file
    file = gzip.open(filename, 'rb')

    basename = filename
    if basename.endswith('.gz'):
        basename = basename[:-3]
    if basename.endswith('.xml'):
        basename = basename[:-4]

    # parse the document
    doc = xml.dom.minidom.parse(file)

    # get the top-level element
    testcase = doc.documentElement

    # open vtf file
    vtffilename = basename + '.vtf'
    vtffile = open(vtffilename, 'w')
    print "Writing configurations to file " + vtffilename + "..."

    configurations = testcase.getElementsByTagName("configuration")
    configuration = configurations[0]
    
    # Output structure
    aid = 0
    for particle in configuration.getElementsByTagName("particle"):
        q = float(particle.attributes['q'].value)
        potential = float(particle.attributes['potential'].value)
        if q < 0.001: 
            name = "N"
        elif q > 0.001:
            name = "O"
        else:
            name = "C"
        vtffile.write("atom {} name {} q {} b {} radius 0.1\n".format(aid, name, q, potential))
        aid = aid + 1

    # Output timesteps
    config_count = 0
    for configuration in configurations:
        print "  %d" % config_count
        
        vtffile.write("\ntimestep\n")

        # unitcell
        box_a = configuration.attributes['box_a'].value.split()
        box_b = configuration.attributes['box_b'].value.split()
        box_c = configuration.attributes['box_c'].value.split()
        vtffile.write("pbc %s %s %s\n" % (box_a[0], box_b[1], box_c[2]))

        # coordinates
        for particle in configuration.getElementsByTagName("particle"):
            vtffile.write(particle.attributes['position'].value + "\n")
            aid = aid + 1

print "Done."
