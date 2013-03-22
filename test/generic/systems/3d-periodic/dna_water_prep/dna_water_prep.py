#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Script that creates a testcase template and includes data from a PDB file.
#
import gzip, sys, re
import xml.dom.minidom

# input files
pdbfilename = 'protonated_imotif_hairpin.pdb'
itpfilename = '1ART6_A.itp'

# output files
vtffilename = 'dna_water.vtf.gz'
xmlfilename = 'dna_water_prep.xml.gz'

template = """<?xml version="1.0" ?>
<!DOCTYPE scafacos_test  SYSTEM 'scafacos_test.dtd'>
<scafacos_test name="dna_water"
               description="A short DNA fragment (i-motif hairpin) in water" 
               reference_method="ewald"
               error_field="1.0e-7" error_potential="1.0e-7" 
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
print('Reading %s...' % pdbfilename)
pdbfile = open(pdbfilename, 'r')
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
print("Read {} atoms from {}.".format(atomcount, pdbfilename))

# READ ITP FILE
qs = []
bonds = []
print("Reading {}...".format(itpfilename))
itpfile = open(itpfilename, 'r')
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
    atoms[no+1]['q'] = 0.417
    atoms[no+2]['q'] = 0.417
    bonds.append((no, no+1))
    bonds.append((no, no+2))

print "Writing to {}...".format(xmlfilename)
xmlfile = gzip.open(xmlfilename, 'w')

sum_q = 0.0
for no in range(1, 15544):
    particle = doc.createElement('particle')
    config_node.appendChild(particle)
    particle.attributes['position'] = "{} {} {}".format(atoms[no]['x'], 
                                                        atoms[no]['y'], 
                                                        atoms[no]['z'])
    particle.attributes['q'] = str(atoms[no]['q'])
    sum_q += atoms[no]['q']

print "sum =", sum_q

xmlfile.write(doc.toprettyxml())
xmlfile.close()

print "Writing to {}...".format(vtffilename)
vtffile = gzip.open(vtffilename, 'w')
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
