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
    vtffilename = basename + '.vtk'
    vtffile = open(vtffilename, 'w')
    print "Writing configurations to file " + vtffilename + "..."

    configurations = testcase.getElementsByTagName("configuration")
    configuration = configurations[0]

    # Output header
    vtffile.write("# vtk DataFile Version 3.0\n")
    vtffile.write("Original file: %s\n" % basename)
    vtffile.write("ASCII\n")
    vtffile.write("DATASET UNSTRUCTURED_GRID\n")

    # Output coordinates
    particles = configuration.getElementsByTagName("particle")
    npart     = len(particles)

    vtffile.write("POINTS %d float\n" % npart)
    for particle in particles:
        vtffile.write(particle.attributes['position'].value + "\n")

    # Output cells and cell types (VTK_VERTEX=1)
#    vtffile.write("CELLS %d %d\n" % (npart, 2*npart))
#    ipart = 0
#    for particle in particles:
#        ipart = ipart + 1
#        vtffile.write("1 %d\n" % ipart)

#    vtffile.write("CELLS_TYPES %d" % npart)
#    for particle in particles:
#        vtffile.write("1\n")

    # Output charges
    vtffile.write("POINT_DATA %d\n" % npart)
    vtffile.write("SCALARS charge float 1\n")
    vtffile.write("LOOKUP_TABLE default\n")
    for particle in particles:
        vtffile.write(particle.attributes['q'].value + "\n")
    
print "Done."
