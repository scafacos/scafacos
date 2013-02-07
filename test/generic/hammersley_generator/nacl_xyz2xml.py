#!/usr/bin/python
#
# Script that creates a testcase template and includes data from an xyz file.
#
import gzip
import sys
import xml.dom.minidom

template = """<?xml version="1.0" ?>
<!DOCTYPE scafacos_test  SYSTEM 'scafacos_test.dtd'>
<scafacos_test name="Name of the testcase"
               description="Short description of the testcase." 
               reference_method="What method was used to generate the reference data."
               error_field="1.0" error_potential="1.0" 
               ></scafacos_test>"""
doc = xml.dom.minidom.parseString(template)

for filename in sys.argv[1:]:
    testcase = doc.documentElement

    config_node = doc.createElement('configuration')
    testcase.appendChild(config_node)

    config_node.attributes['offset'] = '0.0 0.0 0.0'
    config_node.attributes['box_a'] = '0.0 0.0 0.0'
    config_node.attributes['box_b'] = '0.0 0.0 0.0'
    config_node.attributes['box_c'] = '0.0 0.0 0.0'
    config_node.attributes['epsilon'] = 'metallic'
    config_node.attributes['periodicity'] = '1 1 1'

    # open file
    infile = open(filename, 'r')

    madelung = -1.7475645946331821906362120355443974034851614366247417581528253507
    m0 = -1.747564
    m1 = 59463318
    m2 = 21906362
    m3 = 12035544
    m4 = 39740348
    m5 = 51614366
    m6 = 24741758
    m7 = 1528253507
    linecount = 0
    for line in infile.readlines():
        linecount += 1
        line=line.strip()
        
        if len(line) > 0 and not line.startswith('#'):
            try:
                (x,y,z,q) = line.split()
            except ValueError:
                sys.stderr.write("Bad line %d: %s\n" % (linecount, line))
                sys.stderr.write("Should be x y z q\n" % (linecount, line))
            else:
                particle = doc.createElement('particle')
                config_node.appendChild(particle)
                particle.attributes['position'] = "%s %s %s" % (x,y,z)
                particle.attributes['q'] = q
#                particle.attributes['potential'] = "%.17f" % float(q)*madelung
                particle.attributes['potential'] = "%f%d%d%d%d%d%d%d" % (float(q)*m0,m1,m2,m3,m4,m5,m6,m7)
                particle.attributes['field'] = "0 0 0"


print doc.toprettyxml()
