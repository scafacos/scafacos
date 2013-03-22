#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Script that creates a testcase template and includes data from a PDB file.
#
import gzip, sys, re
import xml.dom.minidom

template = """<?xml version="1.0" ?>
<!DOCTYPE scafacos_test  SYSTEM 'scafacos_test.dtd'>
<scafacos_test name="Name of the testcase"
               description="Short description of the testcase." 
               reference_method="What method was used to generate the reference data."
               error_field="1.0" error_potential="1.0" 
               ></scafacos_test>"""
doc = xml.dom.minidom.parseString(template)

testcase = doc.documentElement

config_node = doc.createElement('configuration')
testcase.appendChild(config_node)

config_node.attributes['offset'] = '0.0 0.0 0.0'
config_node.attributes['box_a'] = '54.104 0.0 0.0'
config_node.attributes['box_b'] = '0.0 54.104 0.0'
config_node.attributes['box_c'] = '0.0 0.0 54.104'
config_node.attributes['epsilon'] = 'metallic'
config_node.attributes['periodicity'] = '1 1 1'

# READ PDB FILE    
filename = 'protonated_imotif_hairpin.pdb'
print('Reading %s...' % filename)
pdbfile = open(filename, 'r')
atoms = dict()
linecount = 0
atomcount = 0
for line in pdbfile.readlines():
    linecount += 1
    line=line.strip()

    if line[0:4] == 'ATOM':
        atomcount += 1
        no = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:22].strip()
        matchit = re.match("(.+)\W+(.+)", resname)
        if matchit:
            resname = matchit.group(1)
            chain = matchit.group(2)
        else:
            chain = 'X'
        resid = int(line[22:28])
        x = line[31:38]
        y = line[39:46]
        z = line[47:54]

        atoms[no] = {
            'name': name,
            'resname': resname,
            'resid': resid,
            'chain': chain,
            'x': x,
            'y': y,
            'z': z,
            }

pdbfile.close()
print("Read {} atoms from {}.".format(atomcount, filename))

# READ ITP FILE
filename = '1ART6_A.itp'
qs = []
bonds = []
print("Reading {}...".format(filename))
itpfile = open(filename, 'r')
linecount = 0
while True:
  line = itpfile.readline()
  linecount += 1

  if line == '': break
  line = line.strip()

  if line == '' or line[0] == ';': continue

  # read atoms block
  if line == '[ atoms ]':
        print("  Reading atoms...")
        while True:
            line = itpfile.readline()
            linecount += 1
            if line == '\n': break
            if line[0] == ';': continue
            no = int(line[0:6])
            q = float(line[49:56])
            qs.append(q)
            atoms[no]['q'] = q
  if line == '[ bonds ]':
        print("  Reading bonds...")
        while True:
            line = itpfile.readline()
            linecount += 1
            if line == '\n': break
            if line[0] == ';': continue
            no1 = int(line[0:5])
            no2 = int(line[6:11])
            bonds.append((no1, no2))
itpfile.close()

for no in range(len(qs)):
    atoms[no+1]['q'] = qs[no]

for no in range(687, 703):
    atoms[no]['q'] = 1.0

for no in range(703, 15544, 3):
    atoms[no]['q'] = -0.834
    atoms[no+1]['q'] = -0.417
    atoms[no+2]['q'] = -0.417
    bonds.append((no, no+1))
    bonds.append((no, no+2))

filename = 'dna_water_15543.xml.gz'
print "Writing to {}...".format(filename)
xmlfile = gzip.open(filename, 'w')

for no in range(1, 15544):
    particle = doc.createElement('particle')
    config_node.appendChild(particle)
    particle.attributes['position'] = "{} {} {}".format(atoms[no]['x'], 
                                                        atoms[no]['y'], 
                                                        atoms[no]['z'])
    particle.attributes['q'] = str(atoms[no]['q'])

xmlfile.write(doc.toprettyxml())
xmlfile.close()

filename = 'dna_water_15543.vtf.gz'
print "Writing to {}...".format(filename)
vtffile = gzip.open(filename, 'w')
vtffile.write('pbc 54.104 54.104 54.104\n')
for no in range(1, 15544):
    atom = atoms[no]
    vtffile.write('atom {}'.format(no-1))
    vtffile.write(' name {}'.format(atom['name']))
    vtffile.write(' q {}'.format(atom['q']))
    vtffile.write(' resid {}'.format(atom['resid']))
    vtffile.write(' chain {}'.format(atom['chain']))
    vtffile.write(' resname {}'.format(atom['resname']))
    vtffile.write('\n')
for (p1, p2) in bonds:
    vtffile.write('bond {}:{}\n'.format(p1-1, p2-1))
vtffile.write('timestep\n')
for no in range(1, 15544):
    atom = atoms[no]
    vtffile.write('{} {} {}\n'.format(atom['x'], atom['y'], atom['z']))
vtffile.close()
print "Finished."
