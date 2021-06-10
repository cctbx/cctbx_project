from suitenamedefs import Suite, Residue, globals

"""
This module handles reading suites from "dangle" format files
reading residues from kinemage format files and regrouping them into suites.
Extraction of suites from loaded cctbx models is handled elsewhere.
"""


#                Copyright 2021  Richardson Lab
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import numpy as np
import math, sys


# The great variety of codes that may represent each base in the input file
NAListA = ":ADE:  A:A  : Ar:ATP:ADP:AMP:T6A:1MA:RIA:  I:I  :"
NAListG = ":GUA:  G:G  : Gr:GTP:GDP:GMP:GSP:1MG:2MG:M2G:OMG: 7MG:"
NAListC = ":CYT:  C:C  : Cr:CTP:CDP:CMP:5MC:OMC:"
NAListU = ":URA:URI:  U: Ur:U  :UTP:UDP:UMP:5MU:H2U:PSU:4SU:"
NAListY = ": YG:YG :  Y:Y  :"
#NAListT = ":THY:  T:T  : Tr:TTP:TDP:TMP:"
IgnoreDNAList = ": DA: DG: DC: DT:THY:  T:T  : Tr:TTP:TDP:TMP:"


# out of the noise, determine the base
def findBase(baseCode):
    if len(baseCode) != 3:
        return "Z"

    if NAListA.find(baseCode) >= 0:
        base = "A"
    elif NAListG.find(baseCode) >= 0:
        base = "G"
    elif NAListC.find(baseCode) >= 0:
        base = "C"
    elif NAListU.find(baseCode) >= 0:
        base = "U"
    elif NAListY.find(baseCode) >= 0:
        base = "Y"
    elif IgnoreDNAList.find(baseCode) >= 0:
        return None  # we ignore DNA residues
    else:
        base = "?"
    return base


def stringToFloat(string):
    try:
        n = float(string)
    except ValueError:
        n = 9999.0  # or maybe math.nan?
    return n


def readResidues(inFile):
    options = globals.options

    lines = inFile.readlines()
    residues = []
    for line in lines:
        if len(line.strip()) == 0 or line[0] == "#":  # blank or comment line
            continue
        fields = line.split(":")
        ids = fields[: options.pointidfields]
        # if ids[3].strip() != "":
        #   print('yes')

        baseCode = fields[options.pointidfields - 1]
        angleStrings = fields[options.pointidfields :]
        if (
            ids[options.altidfield].strip() != ""
            and ids[options.altidfield] != options.altidval
        ):
            continue  # lines for the wrong alternative conformation are ignored

        base = findBase(baseCode)
        if not base:  # ignore DNA bases
            continue
        angles = np.array([stringToFloat(s) for s in angleStrings])
        for i in range(len(angles)):
            if angles[i] < 0:
                angles[i] += 360.0

        residue = Residue(ids, base, angles)
        residues.append(residue)
    return residues


def readKinemageFile(inFile):
    """
    We glean the following information from a kinemage file:
    The @dimension command gives us the number of dimensions in the data
    Anything between a @balllist command and a subsequent @ command
    is a data line.
    """
    options = globals.options

    lines = inFile.readlines()
    goodLines = []
    place, line = findPrefixInList(lines, "@dimension")
    if place > 0:
        items = line.split()
        dimension = len(items) - 1
    else:
        dimension = options.anglefields
        place = 0
    while place >= 0:
        begin, line = findPrefixesInList(lines, "@balllist", "@dotlist", place)
        if begin > 0:
            end, line = findPrefixInList(lines, "@", begin + 1)
            place = end
            if end < 0:
                end = len(lines)
            goodLines += lines[begin + 1 : end]
        else:
            break
    if len(goodLines) == 0:
        goodLines = lines  # assume a pure data file
    return readKinemageSuites(goodLines, dimension)


def readKinemageSuites(lines, dimension):
    """Read a list of kinemage data lines to yield a suite."""
    suites = []
    for line in lines:
        if len(line.strip()) == 0 or line[0] == "#":  # blank or comment line
            continue
        # A meaningful line begins with an id string enclosed in braces
        if line[0] == "{":
            mark = line.find("}")
            if mark > 0:
                idString = line[1:mark]
                ids = idString.split(":")

                # there may be some miscellaneous markers after the id string
                k = mark + 1
                while k < len(line) and not line[k].isdigit():
                    k = k + 1
                mark2 = k

                # once we see a number, everything else is angles
                angleText = line[mark2:]
                angleStrings = angleText.split(" ")
                angleStrings2 = angleText.split(",")
                if len(angleStrings2) > len(angleStrings):
                    angleStrings = angleStrings2
                angleList = [stringToFloat(s) for s in angleStrings]
                if len(angleList) != dimension:
                    continue  # wrong number of dimensions means probably not a data point
                if dimension == 9:
                    angles = np.array(angleList)
                else:  # given only 7 angles,skipping the chi angles on the ends
                    angles = np.array([180.0] + angleList + [180.0])
                for i in range(len(angles)):
                    if angles[i] < 0:
                        angles[i] += 360.0

                suite = Suite(ids, ids[9][2], angles)
                suites.append(suite)
    return suites


def findPrefixInList(list, prefix, start=0):
    for i, s in enumerate(list[start:]):
        if s.startswith(prefix):
            return i + start, s
    return -1, None


def findPrefixesInList(list, prefix1, prefix2, start=0):
    for i, s in enumerate(list[start:]):
        if s.startswith(prefix1) or s.startswith(prefix2):
            return i + start, s
    return -1, None


def buildSuiteBetweenResidues(r1, r2):
    suite = Suite(r2.pointIDs, r2.base)
    if len(r1.angle) > 6:
        suite.chiMinus = r1.chi
    suite.deltaMinus = r1.delta
    suite.epsilon = r1.epsilon
    suite.zeta = r1.zeta
    suite.alpha = r2.alpha
    suite.beta = r2.beta
    suite.gamma = r2.gamma
    suite.delta = r2.delta
    if len(r2.angle) > 6:
        suite.chi = r2.chi
    return suite


def buildSuiteFirst(r2):
    suite = Suite(r2.pointIDs, r2.base)
    suite.alpha = r2.alpha
    suite.beta = r2.beta
    suite.gamma = r2.gamma
    suite.delta = r2.delta
    if len(r2.angle) > 6:
        suite.chi = r2.chi
    suite.epsilon = 999
    suite.zeta = 999
    suite.chiMinus = 999
    suite.deltaMinus = 999
    return suite


def buildSuiteLast(r1):
    suite = Suite((), "")
    if len(r1.angle) > 6:
        suite.chiMinus = r1.chi
    suite.deltaMinus = r1.delta
    suite.epsilon = r1.epsilon
    suite.zeta = r1.zeta
    suite.alpha = 999
    suite.beta = 999
    suite.gamma = 999
    suite.delta = 999
    suite.chi = 999
    return suite


def buildSuites(residues):
    suites = [buildSuiteFirst(residues[0])]
    for i in range(len(residues) - 1):
        suites.append(buildSuiteBetweenResidues(residues[i], residues[i + 1]))
    suites.append(buildSuiteLast(residues[-1]))
    return suites
