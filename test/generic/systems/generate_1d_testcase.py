#!/usr/bin/python
# -*- coding: utf-8 -*-
#
from __future__ import print_function
import xml.dom.minidom
import sys
import gzip
from numpy import *

def vanderCorput(N, p):
    # zu wandelnde Zahlen
    numbers = arange(1,int(N)+1)
    # bitumgekehrtes Ergebnis
    result = zeros(N)
    # Wert der aktuellen, inversen Stelle
    frac = 1.0 / p

    # solange die groesste Zahl noch Stellen hat
    while numbers[-1] > 0:
        # unterste Stelle abschneiden
        digit = numbers % p
        numbers /= p
        # ... und zum Ergebnis hinzufuegen
        result += frac*digit
        frac /= p

    return result

try:
    # number of charges
    N = int(sys.argv[1])
except:
    sys.stderr.write('Usage: {} N <filename>\n'.format(sys.argv[0]))
    sys.exit(2)

outname = 'condensator_{}.xml.gz'.format(N)
if len(sys.argv) > 2:
    outname = sys.argv[2]
    if not outname.endswith('.xml.gz'):
        outname += '.xml.gz'

# distance between charges on condensator
dist_cond_charges = 1.0
# minimal distance of ions from condensator
dist_ions = 1.0
# aspect ration between periodic and non-periodic directions
ratio = 1.0
# maximal extent of charges
cloud_oversize = 1.5

# charges per side of condensator
n = sqrt(N/4)
if n.is_integer():
    n = int(n)
else:
    n = int(math.ceil(n))
    n += n%2

# number of ions
Nions = N - 2*n**2

assert(Nions % 2 == 0)
assert(n % 2 == 0)
    
L_xy = n*dist_cond_charges*cloud_oversize
L_z = ratio*L_xy
xoff = 0.5*L_xy - (n-1)*dist_cond_charges/2.0

print('n={} L_xy={} density={}'.format(n, L_xy, N/(L_xy*L_xy*L_z)))

template = """<?xml version="1.0" ?>
<!DOCTYPE scafacos_test  SYSTEM 'scafacos_test.dtd'>
<scafacos_test name="1D periodic testcase"
               description="A periodically replicated plate condensator." 
               reference_method="What method was used to generate the reference data."
               error_field="1.0" error_potential="1.0" 
               ></scafacos_test>"""
doc = xml.dom.minidom.parseString(template)
testcase = doc.documentElement
config_node = doc.createElement('configuration')
testcase.appendChild(config_node)

z_cond1 = L_z*1./4.
q_cond1 = 1.0
z_cond2 = L_z*3./4.
q_cond2 = -1.0

# x/y coordinates of condensator charges
xs = arange(0, n) * dist_cond_charges + xoff
for x in xs:
    for y in xs:
        particle = doc.createElement('particle')
        config_node.appendChild(particle)
        particle.attributes['position'] = "{} {} {}".format(x,y,z_cond1)
        particle.attributes['q'] = str(q_cond1)

for x in xs:
    for y in xs:
        particle = doc.createElement('particle')
        config_node.appendChild(particle)
        particle.attributes['position'] = "{} {} {}".format(x,y,z_cond2)
        particle.attributes['q'] = str(q_cond2)

def create_particles(doc,  config_node,  xs,  ys,  zs,  qs):
    assert(xs.shape == qs.shape == zs.shape == ys.shape)
        
# create ions
xs = vanderCorput(Nions, 3) * L_xy
ys = vanderCorput(Nions, 5) * L_xy
zs = vanderCorput(Nions, 7) * (0.25*L_z - dist_ions)
zs[Nions*1/4:Nions*2/4] += 0.25*L_z + dist_ions
zs[Nions*2/4:Nions*3/4] += 0.5*L_z
zs[Nions*3/4:]          += 0.75*L_z + dist_ions

q = 1.0
for x, y, z in zip(xs, ys, zs):
    particle = doc.createElement('particle')
    config_node.appendChild(particle)
    particle.attributes['position'] = "{} {} {}".format(x,y,z)
    particle.attributes['q'] = str(q)
    q *= -1

config_node.attributes['offset'] = '0.0 0.0 0.0'
config_node.attributes['box_a'] = '{} 0.0 0.0'.format(L_xy)
config_node.attributes['box_b'] = '0.0 {} 0.0'.format(L_xy)
config_node.attributes['box_c'] = '0.0 0.0 {}'.format(L_z)
config_node.attributes['epsilon'] = 'metallic'
config_node.attributes['periodicity'] = '0 0 1'

# Output
print('Writing document to {}...'.format(outname))
outfile = gzip.open(outname, "wb")
outfile.write(doc.toprettyxml())
outfile.close()
print("Finished.")
