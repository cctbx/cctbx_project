# This is an infrequently run file used to generate
# unit test data from the data inherent in suitename.py
# Its output was used to create the canonicalOutput string in UnitTest.py

from __future__ import division
import os, sys, inspect, copy

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from suitename import triageFilters, sieveDelta, sieveGamma
from suiteninit import bins
from suitenamedefs import Residue
from suiteninput import readResidues

import numpy as np
from io import StringIO


def run():
    global pair
    "generate a residue set that will fail each triage in turn"
    s1a = " :1a: : : : Z:   294.967:  173.990:   53.550:   81.495:  212.250:  288.831:  180.000"
    standard = readResidues(StringIO(s1a))
    standard.append(copy.deepcopy(standard[0]))
    # standard is now two copies of the "1a" cluster center
    pair = 2 * standard
    for key, filter in triageFilters.items():
        testFromFilter(key, filter)
    testFromSieve("delta-1", 1, sieveDelta)
    testFromSieve("gamma", 6, sieveGamma)
    testFromSieve("delta", 7, sieveDelta)


def testFromFilter(id, filter):
    index = filter[0]
    min = filter[1]
    testCore(id, index, min)


def testFromSieve(id, index, sieve):
    min = sieve[0][0]
    testCore(id, index, min)


def testCore(id, index, min):
    basic = copy.deepcopy(pair)
    if index >= 4:
        m = index - 4
        n = 1
    else:
        m = index + 2
        n = 0
    basic[n].angle[m] = min - 1.0
    r1 = Residue(" :{}: : : : Z:".format(id), "Z", basic[0].angle)
    r2 = Residue(" :{}: : : : Z:",format(id), "Z", basic[1].angle)
    printResidue(r1)
    printResidue(r2)


def run1():
    "Generate a residue set from the cluster definitions"
    for i in range(1, 13):
        b = bins[i]
        for c in b.cluster[1:]:
            residuesFromCluster(c.angle, c.name)


def residuesFromCluster(angles, id):
    angles1 = np.array(3 * [9999] + list(angles[1:4]) + list(angles[0:1]))
    angles2 = np.array(list(angles[4:8]) + 2 * [9999] + list(angles[8:9]))
    r1 = Residue(" :{}: : : : Z:".format(id), "Z", angles1)
    r2 = Residue(" :{}: : : : Z:".format(id), "Z", angles2)
    printResidue(r1)
    printResidue(r2)


def printResidue(r):
    angleText = ":".join(["{:9.3f}".format(x) for x in r.angle])
    print(r.pointIDs, angleText)


temp = """chiMinus
    # 1   deltaMinus
    # 2   epsilon
    # 3   zeta
    # 4   alpha
    # 5   beta
    # 6   gamma
    # 7   delta
    # 8   chi    """


def rearrange():
    words = []
    for line in temp.splitlines():
        word = line.split()[-1]
        words.append(word)
    text = ", ".join(words)
    return text


run()
