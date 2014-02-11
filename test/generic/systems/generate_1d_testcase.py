#!/usr/bin/python
# -*- coding: utf-8 -*-
#

import xml.dom.minidom
import sys, gzip
from numpy import *

try:
    # number of charges
    N = int(sys.argv[1])
except:
    sys.stderr.write('Usage: {} N <filename>\n'.format(sys.argv[0]))
    sys.exit(2)
    
outname = 'nanowire_{}.xml.gz'.format(N)
if len(sys.argv) > 2:
    outname = sys.argv[2]
    if not outname.endswith('.xml.gz'):
        outname += '.xml.gz'

assert(N%4 == 0)

# line charge density of charges on nano wire
line_density = 0.5
# aspect ratio of the system
aspect_ratio = 1.0
# minimal distance of counterions from wire
r0 = 1.0

template = """<?xml version="1.0" ?>
<!DOCTYPE scafacos_test  SYSTEM 'scafacos_test.dtd'>
<scafacos_test name="1D periodic testcase"
               description="Nanowire." 
               reference_method="What method was used to generate the reference data."
               error_field="1.0" error_potential="1.0" 
               ></scafacos_test>"""
doc = xml.dom.minidom.parseString(template)
testcase = doc.documentElement
config_node = doc.createElement('configuration')
testcase.appendChild(config_node)

z_stride = 1./line_density
L_z = N/2 * z_stride

def create_particles(doc,  config_node,  xs,  ys,  zs,  qs):
    assert(xs.shape == qs.shape == zs.shape == ys.shape)
    for i in range(xs.shape[0]):
        particle = doc.createElement('particle')
        config_node.appendChild(particle)
        particle.attributes['position'] = "{} {} {}".format(xs[i],ys[i],zs[i])
        particle.attributes['q'] = str(qs[i])
        
# create charges on the wire
zs = arange(0, N/2) * z_stride
xs = ys = zeros_like(zs)
qs = ones_like(zs)
qs[N/4:] *= -1
create_particles(doc, config_node, xs, ys, zs, qs)

# create counterions
zs = empty(N/2)
zs[:N/4] = random.uniform(0.0, L_z/2, N/4)
zs[N/4:] = random.uniform(L_z/2, L_z, N/4)
rs = random.exponential(1.0, N/2)
L = aspect_ratio * L_z
rs *= (L-r0)/max(rs)+r0
alphas = random.uniform(0.0, 2*pi, N/2)
xs = rs*sin(alphas)
ys = rs*cos(alphas)
qs *= -1
create_particles(doc, config_node, xs, ys, zs, qs)

# Determine minimal distance between two counterions
# themin = r0
# for i in xrange(len(xs)-1): 
#     dx = xs[i+1:] - xs[i]
#     dy = ys[i+1:] - ys[i]
#     themin = min(themin, min(dx**2+dy**2))
# sys.stderr.write('min_dist={}\n'.format(themin))
    
max_x = max(xs)
min_x = min(xs)
max_y = max(ys)
min_y = min(ys)

L_xy = max(max_x-min_x, max_y-min_y)
sys.stderr.write("L={} {} {}\n".format(L_xy, L_xy, L_z))

config_node.attributes['offset'] = '{} {} 0.0'.format(min_x, min_y)
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
