#!/usr/bin/env python

import sys

class A():
    def __init__(self, l):
        gene, code, change, loc = l.split()
        self.gene = gene
        self.genepos = int(change[:-3])
        self.ref = change[-3]
        self.alt = change[-1]
        self.line = l


with open('/Users/vladsaveliev/vagrant/variantannotation/test/bcra/BRCA_FM.txt') as fm:
    recs = []
    for l in fm:
        recs.append(A(l))

