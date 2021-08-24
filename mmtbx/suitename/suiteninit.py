"""
This is a self initializing module that embodies the data around which
this program is built. It exports its primary data structure:
  bins    the bin and cluster definitions
"""

#        Copyright 2021  Richardson Lab at Duke University
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

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import suitenamedefs
from suitenamedefs import Bin, Cluster, SatelliteInfo

from numpy import array
import argparse, sys

MAX_CLUSTERS = 16  # practical, observed limit of clusters in a bin


# ***parseCommandLine()*******************************************************
def parseCommandLine():
    "Parse a command line, returning a parseargs.Namespace"
    for i, arg in enumerate(sys.argv):
        sys.argv[i] = arg.lower()

    parser = argparse.ArgumentParser()
    buildParser(parser)

    # now actually parse them
    args = parser.parse_args()
    if args.ptid:
        args.pointidfields = args.ptid
    pass
    return args


def buildParser(parser):
    "Set up the expected arguments, ready for parsing"

    # the input file may be given as an argument or as a redirect
    parser.add_argument("infile", nargs="?", default="")

    # input styles
    parser.add_argument("--residuein", "-residuein", action="store_true")
    parser.add_argument("--residuesin", "-residuesin", action="store_true")
    parser.add_argument("--suitein", "-suitein", action="store_true")
    parser.add_argument("--suitesin", "-suitesin", action="store_true")

    # output styles (default is --report)
    outputStyle = parser.add_mutually_exclusive_group()
    outputStyle.add_argument("--report", "-report", action="store_true")
    outputStyle.add_argument("--string", "-string", action="store_true")
    outputStyle.add_argument("--kinemage", "-kinemage", action="store_true")
    # output modifiers
    parser.add_argument("--chart", "-chart", action="store_true")
      # a modifier to --report, suppress as the statistical summary
    parser.add_argument("--causes", "-causes", action="store_true")
      # a modifier to --report, reveals algorithm details
    parser.add_argument("--nosequence", "-nosequence", action="store_true")
      # a modifier to --string, places ':' instead of residue code

    # additional options
    parser.add_argument("--satellites", "-satellites", action="store_true")
    parser.add_argument("--nowannabe", "-nowannabe", action="store_true")
    parser.add_argument("--noinc", "-noinc", action="store_true")
    parser.add_argument("--thetaeta", "-thetaeta", action="store_true")
    parser.add_argument("--etatheta", "-etatheta", action="store_true")
    parser.add_argument("--test", "-test", action="store_true")
    parser.add_argument("--version", "-version", action="store_true")

    # numerical options
    parser.add_argument("--anglefields", "-anglefields", type=int, default=9)
    parser.add_argument("--pointidfields", "-pointidfields", type=int, default=7)
    parser.add_argument("--ptid", "-ptid", type=int, default=0)
    parser.add_argument("--altid", "-altid", type=str, default="A")
    parser.add_argument("--altidval", "-altidval", type=str, default="A")
    parser.add_argument("--altidfield", "-altidfield", type=int, default=6)

    # the following are deprecated:
    parser.add_argument("--angles", type=int, default=9)
    parser.add_argument("--resAngles", type=int, default=6)
    parser.add_argument("--oneline", "-oneline", action="store_true")
    #   "--help" is automatically available, it summarizes this list.
    return parser


# *** codes to match various residues ****************************

idFields = 5
match_list = (
    (":ADE:  A:A  : Ar:ATP:ADP:AMP:T6A:1MA:RIA:  I:I  :", "A"),
    (":GUA:  G:G  : Gr:GTP:GDP:GMP:GSP:1MG:2MG:M2G:OMG: YG: 7MG:YG :", "G"),
    (":CYT:  C:C  : Cr:CTP:CDP:CMP:5MC:OMC:", "C"),
    (":URA:URI:  U: Ur:U  :UTP:UDP:UMP:5MU:H2U:PSU:4SU:", "U"),
    (":THY:  T:T  : Tr:TTP:TDP:TMP:", "T"),
)


# *** Special data relating to satellite clusters ****************

# This function operates on the satelliteData list below
# it creates an associated dictionary based on the name
def buildSatelliteTable():
    global satelliteTable
    satelliteTable = {}
    for item in satelliteData:
        name = item[0]
        satWidths = item[1]
        domWidths = item[2]
        satelliteTable[name] = SatelliteInfo(name, satWidths, domWidths)


# The satellite data:
# The widths below are used to determine multidimensional hyperellipsoidal distances.
# A distance <= 1  is considered "in" the cluster

# There are three tiers of widths here:
# 1. The normal widths, deltamw etc.
# 2. The general satellite widths, epsilonsatw etc.
# 3. The special satellite widths in the table farther below, satelliteData
# The normalWidths are used for a typical cluster.
# The satelliteWidths are an exception, used for satellite clusters.
# The satelliteData are exceptions to the satelliteWidths, for certain specific
# satellite clusters.

# SITUATION as of 210213:
#   The satelliteWidths are used ONLY if the --satellites arg is used

clusterhalfwidthsversion = "070328"
deltamw = 28
epsilonw = 60
epsilonsatw = 50  # satw 070328
zetaw = 55
zetasatw = 50  # satw 070328
alphaw = 50
alphasatw = 45  # satw 070328
betaw = 70
betasatw = 60  # satw 070328
gammaw = 35
deltaw = 28


# width arrays set the widths of clusters in the various dimensions
# the zeroes on either end may someday be replaced with widths for the
# chi angles
normalWidths = array((0, deltamw, epsilonw, zetaw, alphaw, betaw, gammaw, deltaw, 0))
satelliteWidths = array(
    (0, deltamw, epsilonsatw, zetasatw, alphasatw, betasatw, gammaw, deltaw, 0)
)

satelliteData = (
    #  sat         9 angles sat widths                    9 angles dom width
    ("1m", (0, 0, 0, 0, 0, 32, 0, 0, 0), (0, 0, 0, 0, 0, 64, 0, 0, 0)),
    ("1L", (0, 0, 18, 0, 0, 18, 0, 0, 0), (0, 0, 70, 0, 0, 70, 0, 0, 0)),
    ("&a", (0, 0, 20, 20, 0, 0, 0, 0, 0), (0, 0, 60, 60, 0, 0, 0, 0, 0)),
    ("1f", (0, 0, 0, 0, 0, 47, 0, 0, 0), (0, 0, 0, 0, 0, 65, 0, 0, 0)),
    ("1[", (0, 0, 0, 0, 0, 34, 0, 0, 0), (0, 0, 0, 0, 0, 56, 0, 0, 0)),
    ("4a", (0, 0, 40, 40, 0, 0, 0, 0, 0), (0, 0, 50, 50, 0, 0, 0, 0, 0)),
    ("#a", (0, 0, 26, 26, 0, 0, 0, 0, 0), (0, 0, 36, 36, 0, 0, 0, 0, 0)),
    ("0i", (0, 0, 0, 0, 0, 60, 0, 0, 0), (0, 0, 0, 0, 0, 60, 0, 0, 0)),
    ("6j", (0, 0, 0, 0, 0, 60, 0, 0, 0), (0, 0, 0, 0, 0, 60, 0, 0, 0)),
)
# note on the satelliteData:
# The satellite cluster is stated. The dominant cluster can be found by looking
# in the bin containing the satellite cluster, it will be the cluster with a
# dominance of "dom"


def getSatelliteInfo(name):
  if name in satelliteTable:
    return satelliteTable[name]
  else:
    return None


# *** Cluster data  *******************************************************

# The cluster data: centers of each cluster in 7 dimensions
#   (number, name, status, color, dominance ... the 7 angles)
bin0data = (0, "trig",
    ( 0 , "!!", "triaged", "white      ", "tri",
        (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
)

bin1data = (1, "33 p",
    ( 0 , "!!", "outlier", "white      ", "out",
        (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
    ( 1 , "1a", "certain", "yellowtint ", "dom",
        (180.0,  81.495,  212.25,  288.831,  294.967,  173.99,  53.55,  81.035,  180.0)),
    ( 2 , "1m", "certain", "blue       ", "sat",
        (180.0,  83.513,  218.12,  291.593,  292.247,  222.3,  58.067,  86.093,  180.0)),
    ( 3 , "1L", "certain", "green      ", "sat",
        (180.0,  85.664,  245.014,  268.257,  303.879,  138.164,  61.95,  79.457,  180.0)),
    ( 4 , "&a", "certain", "cyan       ", "sat",
        (180.0,  82.112,  190.682,  264.945,  295.967,  181.839,  51.455,  81.512,  180.0)),
    ( 5 , "7a", "certain", "pink       ", "ord",
        (180.0,  83.414,  217.4,  222.006,  302.856,  160.719,  49.097,  82.444,  180.0)),
    ( 6 , "3a", "certain", "magenta    ", "ord",
        (180.0,  85.072,  216.324,  173.276,  289.32,  164.132,  45.876,  84.956,  180.0)),
    ( 7 , "9a", "certain", "hotpink    ", "ord",
        (180.0,  83.179,  210.347,  121.474,  288.568,  157.268,  49.347,  81.047,  180.0)),
    ( 8 , "1g", "certain", "sea        ", "ord",
        (180.0,  80.888,  218.636,  290.735,  167.447,  159.565,  51.326,  85.213,  180.0)),
    ( 9 , "7d", "certain", "purple     ", "ord",
        (180.0,  83.856,  238.75,  256.875,  69.562,  170.2,  52.8,  85.287,  180.0)),
    ( 10 , "3d", "certain", "peach      ", "ord",
        (180.0,  85.295,  244.085,  203.815,  65.88,  181.13,  54.68,  86.035,  180.0)),
    ( 11 , "5d", "certain", "yellow     ", "ord",
        (180.0,  79.671,  202.471,  63.064,  68.164,  143.45,  49.664,  82.757,  180.0)),
    ( 12 , "3g", "wannabe", "gray       ", "ord",
        (180.0,  84.0,  195.0,  146.0,  170.0,  170.0,  52.0,  84.0,  180.0)),
)

bin2data = (2, "33 t",
    ( 0 , "!!", "outlier", "white      ", "out",
        (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
    ( 1 , "1e", "certain", "red        ", "ord",
        (180.0,  80.514,  200.545,  280.51,  249.314,  82.662,  167.89,  85.507,  180.0)),
    ( 2 , "1c", "certain", "gold       ", "dom",
        (180.0,  80.223,  196.591,  291.299,  153.06,  194.379,  179.061,  83.648,  180.0)),
    ( 3 , "1f", "certain", "lime       ", "sat",
        (180.0,  81.395,  203.03,  294.445,  172.195,  138.54,  175.565,  84.47,  180.0)),
    ( 4 , "5j", "certain", "sky        ", "ord",
        (180.0,  87.417,  223.558,  80.175,  66.667,  109.15,  176.475,  83.833,  180.0)),
    ( 5 , "5n", "wannabe", "gray       ", "ord",
        (180.0,  86.055,  246.502,  100.392,  73.595,  213.752,  183.395,  85.483,  180.0)),
)

bin3data = (3, "33 m",
    ( 0 , "!!", "outlier", "white      ", "out",
        (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
    # ( 1 , "!!", "nothing", "white      ", "out",
    #     (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0)),
    # viewed as highly dubious, KPB 210308
)

bin4data = (4, "32 p",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "1b", "certain", "cyan       ", "dom",
        (180.000, 084.215, 215.014, 288.672, 300.420, 177.476, 058.307, 144.841, 180.000)),
    ( 2 , "1[", "certain", "pink       ", "sat",
        (180.000, 082.731, 220.463, 288.665, 296.983, 221.654, 054.213, 143.771, 180.000)),
    ( 3 , "3b", "certain", "lilac      ", "ord",
        (180.000, 084.700, 226.400, 168.336, 292.771, 177.629, 048.629, 147.950, 180.000)),
    ( 4 , "1z", "certain", "peach      ", "ord",
        (180.000, 083.358, 206.042, 277.567, 195.700, 161.600, 050.750, 145.258, 180.000)),
    ( 5 , "5z", "certain", "purple     ", "ord",
        (180.000, 082.614, 206.440, 052.524, 163.669, 148.421, 050.176, 147.590, 180.000)),
    ( 6 , "7p", "certain", "sea        ", "ord",
        (180.000, 084.285, 236.600, 220.400, 068.300, 200.122, 053.693, 145.730, 180.000)),
    ( 7 , "5p", "wannabe", "gray       ", "ord",
        (180.000, 084.457, 213.286, 069.086, 075.500, 156.671, 057.486, 147.686, 180.000)),
)

bin5data = (5, "32 t",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "1t", "certain", "red        ", "ord",
        (180.000, 081.200, 199.243, 288.986, 180.286, 194.743, 178.200, 147.386, 180.000)),
    ( 2 , "5q", "certain", "yellow     ", "ord",
        (180.000, 082.133, 204.933, 069.483, 063.417, 115.233, 176.283, 145.733, 180.000)),
)

bin6data = (6, "32 m",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "1o", "certain", "sky        ", "ord",
        (180.000, 083.977, 216.508, 287.192, 297.254, 225.154, 293.738, 150.677, 180.000)),
    ( 2 , "7r", "certain", "lilactint  ", "ord",
        (180.000, 084.606, 232.856, 248.125, 063.269, 181.975, 295.744, 149.744, 180.000)),
    ( 3 , "5r", "wannabe", "gray       ", "ord",
        (180.000, 083.000, 196.900, 065.350, 060.150, 138.425, 292.550, 154.275, 180.000)),
)

bin7data = (7, "23 p",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "2a", "certain", "cyan       ", "ord",
        (180.000, 145.399, 260.339, 288.756, 288.444, 192.733, 053.097, 084.067, 180.000)),
    ( 2 , "4a", "certain", "yellow     ", "sat",
        (180.000, 146.275, 259.783, 169.958, 298.450, 169.583, 050.908, 083.967, 180.000)),
    ( 3 , "0a", "certain", "green      ", "dom",
        (180.000, 149.286, 223.159, 139.421, 284.559, 158.107, 047.900, 084.424, 180.000)),
    ( 4 , "#a", "certain", "hotpink    ", "sat",
        (180.000, 148.006, 191.944, 146.231, 289.288, 150.781, 042.419, 084.956, 180.000)),
    ( 5 , "4g", "certain", "greentint  ", "ord",
        (180.000, 148.028, 256.922, 165.194, 204.961, 165.194, 049.383, 082.983, 180.000)),
    ( 6 , "6g", "certain", "gold       ", "ord",
        (180.000, 145.337, 262.869, 079.588, 203.863, 189.688, 058.000, 084.900, 180.000)),
    ( 7 , "8d", "certain", "red        ", "ord",
        (180.000, 148.992, 270.596, 240.892, 062.225, 176.271, 053.600, 087.262, 180.000)),
    ( 8 , "4d", "certain", "sky        ", "ord",
        (180.000, 149.822, 249.956, 187.678, 080.433, 198.133, 061.000, 089.378, 180.000)),
    ( 9 , "6d", "certain", "orange     ", "ord",
        (180.000, 146.922, 241.222, 088.894, 059.344, 160.683, 052.333, 083.417, 180.000)),
    ( 10 , "2g", "wannabe", "gray       ", "ord",
        (180.000, 141.900, 258.383, 286.517, 178.267, 165.217, 048.350, 084.783, 180.000)),
)

bin8data = (8, "23 t",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "2h", "certain", "sea        ", "ord",
        (180.000, 147.782, 260.712, 290.424, 296.200, 177.282, 175.594, 086.565, 180.000)),
    ( 2 , "4n", "certain", "peach      ", "ord",
        (180.000, 143.722, 227.256, 203.789, 073.856, 216.733, 194.444, 080.911, 180.000)),
    ( 3 , "0i", "certain", "lilactint  ", "sat",
        (180.000, 148.717, 274.683, 100.283, 080.600, 248.133, 181.817, 082.600, 180.000)),
    ( 4 , "6n", "certain", "lilac      ", "dom",
        (180.000, 150.311, 268.383, 084.972, 063.811, 191.483, 176.644, 085.600, 180.000)),
    ( 5 , "6j", "certain", "purple     ", "sat",
        (180.000, 141.633, 244.100, 066.056, 071.667, 122.167, 182.200, 083.622, 180.000)),
)

bin9data = (9, "23 m",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "0k", "wannabe", "gray       ", "ord",
        (180.000, 149.070, 249.780, 111.520, 278.370, 207.780, 287.820, 086.650, 180.000)),
)

bin10data = (10, "22 p",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "2[", "certain", "sea        ", "ord",
        (180.000, 146.383, 259.402, 291.275, 291.982, 210.048, 054.412, 147.760, 180.000)),
    ( 2 , "4b", "certain", "gold       ", "ord",
        (180.000, 145.256, 244.622, 162.822, 294.159, 171.630, 045.900, 145.804, 180.000)),
    ( 3 , "0b", "certain", "red        ", "ord",
        (180.000, 147.593, 248.421, 112.086, 274.943, 164.764, 056.843, 146.264, 180.000)),
    ( 4 , "4p", "certain", "purple     ", "ord",
        (180.000, 150.077, 260.246, 213.785, 071.900, 207.638, 056.715, 148.131, 180.000)),
    ( 5 , "6p", "certain", "sky        ", "ord",
        (180.000, 146.415, 257.831, 089.597, 067.923, 173.051, 055.513, 147.623, 180.000)),
    ( 6 , "2z", "wannabe", "gray       ", "ord",
        (180.000, 142.900, 236.550, 268.800, 180.783, 185.133, 054.467, 143.350, 180.000)),
)

bin11data = (11, "22 t",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "4s", "certain", "lime       ", "ord",
        (180.000, 149.863, 247.562, 170.488, 277.938, 084.425, 176.413, 148.087, 180.000)),
    ( 2 , "2u", "wannabe", "gray       ", "ord",
        (180.000, 143.940, 258.200, 298.240, 279.640, 183.680, 183.080, 145.120, 180.000)),
)

bin12data = (12, "22 m",
    ( 0 , "!!", "outlier", "white      ", "out",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
    ( 1 , "2o", "certain", "hotpink    ", "ord",
        (180.000, 147.342, 256.475, 295.508, 287.408, 194.525, 293.725, 150.458, 180.000)),
)

bin13data = (13, "inc ",
    ( 0 , "__", "incompl", "white      ", "inc",
        (000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000, 000.000)),
)


def buildBin(data):
    ordinal = data[0]
    name = data[1]
    clusters = []
    for item in data[2:]:
        c = suitenamedefs.Cluster(*item)
        if c.dominance == "sat":
            c.satelliteInfo = getSatelliteInfo(c.name)
        clusters.append(c)
    bin = Bin(ordinal, name, clusters)
    return bin


# The bins become an associative table that maps the bin selector information
# directly to the bin
# selector = (puckerdm, puckerd, gammaname)
# bins 0 and 13 are catchbasins for outliers, they are not indexed by selectors
def buildTheBins():
    bins = {}
    bins[0] = buildBin(bin0data)
    bins[(3, 3, "p")] = buildBin(bin1data)
    bins[(3, 3, "t")] = buildBin(bin2data)
    bins[(3, 3, "m")] = buildBin(bin3data)
    bins[(3, 2, "p")] = buildBin(bin4data)
    bins[(3, 2, "t")] = buildBin(bin5data)
    bins[(3, 2, "m")] = buildBin(bin6data)
    bins[(2, 3, "p")] = buildBin(bin7data)
    bins[(2, 3, "t")] = buildBin(bin8data)
    bins[(2, 3, "m")] = buildBin(bin9data)
    bins[(2, 2, "p")] = buildBin(bin10data)
    bins[(2, 2, "t")] = buildBin(bin11data)
    bins[(2, 2, "m")] = buildBin(bin12data)
    bins[13] = buildBin(bin13data)

    # build aliases so that bins can be indexed by number during output
    bins[1] = bins[(3, 3, "p")]
    bins[2] = bins[(3, 3, "t")]
    bins[3] = bins[(3, 3, "m")]
    bins[4] = bins[(3, 2, "p")]
    bins[5] = bins[(3, 2, "t")]
    bins[6] = bins[(3, 2, "m")]
    bins[7] = bins[(2, 3, "p")]
    bins[8] = bins[(2, 3, "t")]
    bins[9] = bins[(2, 3, "m")]
    bins[10] = bins[(2, 2, "p")]
    bins[11] = bins[(2, 2, "t")]
    bins[12] = bins[(2, 2, "m")]
    return bins


# parseCommandLine()
buildSatelliteTable()
bins = buildTheBins()
