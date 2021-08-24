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
from __future__ import  with_statement, print_function, unicode_literals
from suitenamedefs import globals, reasons
from suiteninit import bins

import numpy as np

options = None  # global,  will be copied from suitename.py global

reportCountAll = 0
trigCountAll = 0
suitenessSumAll = 0
binSuiteCountAll = 0


def output(outFile, suites, outNote):
  global options
  options = globals.options
  for s in suites:
    outSuite(outFile, s)
  writeFinalOutput(outFile, suites, outNote)


def outSuite(outFile, s):
  if options.string:
    string1Suite(outFile, s)
  elif options.kinemage:
    pass
  else:
    reportSuite(outFile, s)


def writeFinalOutput(outFile, suites, outNote):
  if options.satellites:
    outNote.comment = " special general case satellite widths, power = 3.00"
  else:
    outNote.comment = " all general case widths, power = 3.00"
  if options.string:
    outFile.write("\n")
  elif options.kinemage:
    kinemageFinal(outFile, suites, outNote)
  else:
    reportFinal(outFile, outNote)


def string1Suite(outFile, suite):
  if options.nosequence:
    basestring = ":"
  else:
    basestring = suite.base
  if suite.cluster:
    name = suite.cluster.name
  else:
    name = "!!"
  outFile.write("{}{}".format(name,basestring))


def reportSuite(outFile, suite):
  global reportCountAll, trigCountAll, suitenessSumAll, binSuiteCountAll

  if options.noinc and not suite.valid:
    return

  # 1. write one line of output for this suite
  reason = ""; note=""
  outIDs = ":".join(suite.pointID)
  if suite.issue:
    reason = " " + reasons[suite.issue]
  elif suite.comment:
    reason = " " + suite.comment
  elif suite.situation and options.causes:
    reason += " " + suite.situation
  if suite.cluster.status == "wannabe":
    note = " wannabe"
  output = (
    "{} {} {} {:5.3f}{}{}\n".format(outIDs, suite.bin.name, suite.cluster.name, float(suite.suiteness), reason, note)
  )
  outFile.write(output)

  # 2. gather statistics
  reportCountAll += 1
  bin = suite.bin
  cluster = suite.cluster
  suiteness = suite.suiteness
  if bin.ordinal == 0:
    trigCountAll += 1
  elif bin.ordinal < 13:
    suitenessSumAll += suiteness
    binSuiteCountAll += 1

  if cluster.ordinal == 0:
    cluster.suitenessCounts[11] += 1
  else:
    cluster.suitenessSum += suiteness
    # report in statistical baskets at intervals of 0.1:
    # everything from 0 to 0.1 goes in bucket 1
    # ... everything from 0.9 to 1.o goes into bucket 10
    if suiteness == 0:
      bucket = 0
    else:
      bucket = 1 + int(suiteness * 10)
    cluster.suitenessCounts[bucket] += 1


def reportFinal(outFile, outNote):
  if not options.chart:
    outFile.write(outNote.comment + "\n")
    suitenessAverage(outFile, 0)
    if bins[1].cluster[1].count > 0:  # Aform 1a    070325
      suitenessAverage(outFile, 1)
      suitenessAverage(outFile, 2)


def suitenessAverage(outFile,  mode):
  """
  Gather statistics on suiteness.
  12 buckets:
    one for suiteness=0
    one each for divisions by tenths from 0 to 1
    one for outliers
  """
  bucket = np.zeros(12, dtype=int)
  sum = 0
  average = 0
  allCount = 0
  excludedCluster = None

  if mode == 1:
    # cluster 1a all by itself
    comment = " A form (1a)"
    cluster = bins[1].cluster[1]
    sum = cluster.suitenessSum
    for k in range(12):
      bucket[k] = cluster.suitenessCounts[k]
  else:
    if mode == 0:
      # all complete suites
      startWith = 0
      comment = "For all"
      outFile.write(
          "Found {} complete suites derived from {} entries\n".format(
          binSuiteCountAll + trigCountAll, reportCountAll
        )
      )
      outFile.write(
          "{} suites were triaged, leaving {} assigned to bins\n".format(
          trigCountAll, binSuiteCountAll
        )
      )
    elif mode == 2:
      # all complete suites except cluster 1a
      startWith = 1
      comment = " non-1a  has"
      excludedCluster = bins[1].cluster[1]

    for i in range(1, 13):
      bin = bins[i]
      # the final bin 13 is for the pseudo-suites with incomplete
      # angles and is ignored
      for cluster in bin.cluster[startWith:]:
        # cluster 0 in every bin is for outliers and is ignored
        if cluster is excludedCluster:
          continue  # mode 2: ignore cluster 1a
        sum += cluster.suitenessSum
        for k in range(12):
          bucket[k] += cluster.suitenessCounts[k]

  allCount = np.sum(bucket)
  if allCount > 1:
    average = sum / allCount
  else:
    average = 0
  outFile.write(
    "{} {} suites: average suiteness== {:5.3f} (power==3.00)\n".format(
        comment, allCount, average
    )
  )
  if mode == 0:
    outFile.write("{:6d} suites are  outliers\n".format(bucket[11]))
  outFile.write("{:6d} suites have suiteness == 0    \n".format(bucket[0]))
  outFile.write("{:6d} suites have suiteness >  0 <.1\n".format(bucket[1]))
  for k in range(2, 10):
    outFile.write(
      "{:6d} suites have suiteness >=.{} <.{}\n".format(bucket[k], k - 1, k)
    )
  outFile.write("{:6d} suites have suiteness >=.9    \n".format(bucket[10]))


def clearStatistics():
  global reportCountAll, trigCountAll, suitenessSumAll, binSuiteCountAll
  reportCountAll = 0
  trigCountAll = 0
  suitenessSumAll = 0
  binSuiteCountAll = 0

  for bin in bins.values():
    for cluster in bin.cluster:
      cluster.count = 0
      cluster.suitenessSum = 0
      for k in range(12):
        cluster.suitenessCounts[k] = 0


# ***** kinemage output format *****************************************************

def kinemage1Suite(suite, bin, cluster, notes, distance, suiteness, issue, comment,
                    pointMaster, pointColor):
    suite.pointMaster = pointMaster
    suite.pointColor = pointColor
    suite.notes = notes


# static text: viewing parameters
janesviews = """
@viewid {d e z}
@zoom 1.00
@zslab 200
@ztran 0
@center 197.500 172.300 178.300
@axischoice 2 3 4
@matrix
0.07196 0.11701 -0.99052 -0.00336 0.99312 0.11707 0.99740 -0.00509 0.07186
@2viewid {zag front}
@2zoom 1.00
@2zslab 200
@2ztran 0
@2center 174.091 194.887 207.768
@2axischoice 4 5 7
@2matrix
0.99508 -0.00018 -0.09905 -0.00135 -0.99993 -0.01172 -0.09904 0.0118 -0.99501
@3viewid {a b g}
@3zoom 1.00
@3zslab 200
@3ztran 0
@3center 175.700 189.600 64.100
@3axischoice 5 6 7
@3matrix
0.99955 0.000101 0.030002 0.0002 0.99995 -0.010012 -0.030001 0.010013 0.9995

"""

# static text: static items in the display
kinemageFrame = ("\n"+
  "@group {frame} dominant \n"
  "@vectorlist {frame} color= white \n"
  "P   0.000   0.000   0.000   5.000   0.000   0.000 \n"+
  "P  35.000   0.000   0.000  40.000   0.000   0.000 \n"+
  "P  80.000   0.000   0.000 160.000   0.000   0.000 \n"+
  "P   0.000   0.000   0.000   0.000   5.000   0.000 \n"+
  "P   0.000  35.000   0.000   0.000  40.000   0.000 \n"+
  "P   0.000  80.000   0.000   0.000 160.000   0.000 \n"+
  "P   0.000   0.000   0.000   0.000   0.000   5.000 \n"+
  "P   0.000   0.000  35.000   0.000   0.000  40.000 \n"+
  "P   0.000   0.000  80.000   0.000   0.000 160.000 \n"+
  "P 200.000   0.000   0.000 280.000   0.000   0.000 \n"+
  "P 320.000   0.000   0.000 360.000   0.000   0.000 \n"+
  "P   0.000 200.000   0.000   0.000 280.000   0.000 \n"+
  "P   0.000 320.000   0.000   0.000 360.000   0.000 \n"+
  "P   0.000   0.000 200.000   0.000   0.000 280.000 \n"+
  "P   0.000   0.000 320.000   0.000   0.000 360.000 \n"+
  "@labellist {XYZ} color= white \n"+
  "{X}  20.000  -5.000  -5.000 \n"+
  "{X} 380.000  -5.000  -5.000 \n"+
  "{Y}  -5.000  20.000  -5.000 \n"+
  "{Y}  -5.000 380.000  -5.000 \n"+
  "{Z}  -5.000  -5.000  20.000 \n"+
  "{Z}  -5.000  -5.000 380.000 \n"+
  "@labellist {mtp} color= green \n"+
  "{p}  60.000   0.000   0.000 \n"+
  "{t} 180.000   0.000   0.000 \n"+
  "{m} 300.000   0.000   0.000 \n"+
  "{p}   0.000  60.000   0.000 \n"+
  "{t}   0.000 180.000   0.000 \n"+
  "{m}   0.000 300.000   0.000 \n"+
  "{p}   0.000   0.000  60.000 \n"+
  "{t}   0.000   0.000 180.000 \n"+
  "{m}   0.000   0.000 300.000 \n"+
  "\n"
)

def kinemageFinal(outFile, suites, outNote):
  """
  Output the content of a kinemage file
  The 3, 2 order may seem odd, but it is a standard
  """
  if not suites or all([not s.valid for s in suites]):
    # empty input gives empty output
    return
  kinemageHeader(outFile, outNote)
  for deltaMinus in (3, 2):
    for delta in (3, 2):
      binGroupOut(outFile, deltaMinus, delta, suites)
  triaged = bins[0]
  if triaged.count > 0:
    outFile.write("@group {triaged} dominant dimension=9 wrap=360 select off\n")
    binOut(outFile, triaged, suites)


def kinemageHeader(outFile, outNote):
  """ The invariant portion of a kinemage file """
  outFile.write("@text\n {}\n {}\n".format(outNote.version, outNote.comment))
  outFile.write("@kinemage 1\n")
  outFile.write("@onewidth\n")
  if options.etatheta: # 070524
    outFile.write(
        "@dimension {theta} {delta-1} {epsilon-1} {zeta-1} {alpha} "
        "{beta} {gamma} {delta} {eta}\n"
    )
  else:
    outFile.write(
        "@dimension {chi-1} {delta-1} {epsilon-1} {zeta-1} {alpha} "
        "{beta} {gamma} {delta} {chi}\n"
    )
  outFile.write(
      "@dimminmax 0.000 360.000 0.000 360.000 0.000 360.000 0.000 "
      "360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 "
      "0.000 360.000\n"
    )
  if outNote.outliers:
    outFile.write("@pointmaster 'O' {outliers}\n")
  if outNote.wannabes:
    outFile.write("@master {wannabees}\n")
  outFile.write(janesviews)
  outFile.write(kinemageFrame)


def binGroupOut(outFile, deltaMinus, delta, suites):
  # If any bin in the group has data, generate a group header
  groupCount = 0
  for gamma in ("p", "t", "m"):
    bin = bins[(deltaMinus, delta, gamma)]
    groupCount += bin.count
  if groupCount > 0:
    outFile.write(
      "@group {{{}{}}} recessiveon dimension=9".format(deltaMinus,delta)+
      " wrap=360 select animate off\n")

  # generate the data
  for gamma in ("p", "t", "m"):
    bin = bins[(deltaMinus, delta, gamma)]
    if bin.count > 0:
      binOut(outFile, bin, suites)


def binOut(outFile, bin, suites):
  if any([c.count > 0 for c in bin.cluster]):
    outFile.write("@subgroup {{{}}} recessiveon \n".format(bin.name))
  for cluster in bin.cluster[1:]:
      # the first cluster, for outliers, will be handled later
      if cluster.count > 0:  # empty clusters vanish
        extras = ""
        if cluster.status == "wannabe":
          extras = " master= {wannabees}"
        # display a ball for each for each suite in this cluster
        ballList = (
            "@balllist {{{} {}}} color= {} radius= 1 "
            "nohilite master= {{data}}{}\n"
        ).format(bin.name, cluster.name, cluster.clusterColor, extras)
        # display a ring surrounding the center of the cluster
        ringList = (
            "@ringlist {{{} {}}} color= {} radius= 10 width= 1 "
            "nobutton master= {{avsigma}}{}\n"
        ).format(bin.name, cluster.name, cluster.clusterColor, extras)
        angleList = formatAngles(cluster.angle[1:-1], ' ')
        ringList2 = "{{{} {}}} 180 {} 180\n".format(
            bin.name, cluster.name, angleList
        )
        labelList = (
            "@labellist {{{} {}}} color= {} nobutton "
            "master= {{labels}}{}\n"
        ).format(bin.name, cluster.name, cluster.clusterColor, extras)
        outFile.write(ballList)
        outPoints(outFile, bin, cluster, suites, "")
        outFile.write(ringList)
        outFile.write(ringList2)
        outFile.write(labelList)
        outFile.write(ringList2)

  # handle outliers if there are any:
  if bin.cluster[0].count > 0:
    cluster = bin.cluster[0]
    ballList = (
      "@balllist {{{} {}}} color= {} radius= 1 "
      "nohilite master= {{data}}\n"
    ).format(bin.name, cluster.name, cluster.clusterColor)
    outFile.write(ballList)
    outPoints(outFile, bin, cluster, suites, "'O' white")


def outPoints(outFile, bin, cluster, suites, extra1):
  for s in suites:
    if s.cluster is cluster:
      if bin.name == "trig":  # the triage bin is specially handled
        extra = "'{}'".format(s.pointMaster)
      else:
        extra = extra1
      ids = ':'.join(s.pointID)
      line = \
        ("{{{} {} :D=={:5.3f}".format(bin.name, cluster.name, s.distance) \
       + ":S=={:5.3f}: {}}} {} ,".format(s.suiteness, ids, extra)) \
       + formatAngles(s.angle, ",") + "\n"
      outFile.write(line)


def formatAngles(angles, separator):
    strings = ["{:7.2f}".format(a) for a in angles]
    out = separator.join(strings)
    return out
