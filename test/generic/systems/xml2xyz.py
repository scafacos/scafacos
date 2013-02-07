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
    vtffilename = basename + '.xyz'
    vtffile = open(vtffilename, 'w')
    print "Writing configurations to file " + vtffilename + "..."

    configurations = testcase.getElementsByTagName("configuration")
    configuration = configurations[0]
    
    # Output structure
    aid = 0
    for particle in configuration.getElementsByTagName("particle"):
        vtffile.write(particle.attributes['position'].value + " " + particle.attributes['q'].value + "\n")
        
print "Done. Visualize your data with gnuplot> splot " + vtffilename + " using 1:2:3:4"
