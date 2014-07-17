import sys
import xml.etree.ElementTree as ElementTree
import gzip

if len(sys.argv) != 4:
    sys.stderr.write("usage: {} <input> \"<srcx> <srcy> <srcz>\" <output>\n".format(sys.argv[0]))
    exit(-1)

shuffle_pattern = [int(c) for c in sys.argv[2].split()]

def shuffle(v):
    v = v.split()
    return " ".join([v[shuffle_pattern[i]] for i in xrange(3)]) 
    
tree = ElementTree.parse(gzip.open(sys.argv[1]))
root = tree.getroot()

for element in root.iter("configuration"):
    element.attrib["periodicity"] = shuffle(element.attrib["periodicity"])
    element.attrib["offset"] = shuffle(element.attrib["offset"])
    # adapt both the coordinates and the order of the lattice vectors, so that uses_principal_axes does not complain
    box = shuffle(element.attrib["box_a"]), shuffle(element.attrib["box_b"]), shuffle(element.attrib["box_c"])
    box = [box[shuffle_pattern[i]] for i in xrange(3)]
    element.attrib["box_a"], element.attrib["box_b"], element.attrib["box_c"] = box
    
for element in root.iter("particle"):
    element.attrib["position"] = shuffle(element.attrib["position"])
    element.attrib["field"] = shuffle(element.attrib["field"])
    
tree.write(sys.argv[3])
