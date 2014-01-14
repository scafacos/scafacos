#!/usr/bin/python
# -*- coding: utf-8 -*-
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
    sys.stderr.write('Reading %s...\n' % filename)
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

    linecount = 0
    for line in infile.readlines():
        linecount += 1
        line=line.strip()
        
        if len(line) > 0 and not line.startswith('#'):
            try:
                (x,y,z,q,p,fx,fy,fz) = line.split()
            except ValueError:
                sys.stderr.write("Bad line %d: %s\n" % (linecount, line))
                sys.stderr.write("Should be x y z q p fx fy fz\n" % (linecount, line))
            else:
                particle = doc.createElement('particle')
                config_node.appendChild(particle)
                particle.attributes['position'] = "%s %s %s" % (x,y,z)
                particle.attributes['q'] = q
                particle.attributes['xpotential'] = p
                particle.attributes['zfield'] = "%s %s %s" % (fx,fy,fz)

print doc.toprettyxml()
