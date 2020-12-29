from xml.dom import minidom
import xml.etree.ElementTree as ET
import numpy as np


with open('test.mzML', 'rt') as f:
   tree = ET.parse(f)

nt = tree.getiterator(tag = 'spectrum')
nt
list(tree.iter())

for node in tree.getroot()
   print (node.tag)

tree.findall('indexedmzML')
