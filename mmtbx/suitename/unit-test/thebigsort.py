# thebigsort.py: Collect data on all the reference database,
# sorted by "situation"

import os, sys
print(os.getcwd())
print(sys.path)
sys.path.insert(0, '.') 
# do this before importing anything from suite name:
sys.argv.extend(("--pointidfields", "7"))

from pathlib import Path
import suitename
from suitenamedefs import reasons

def setup():
    workDir = os.path.dirname(sys.argv[0]).strip()
    if workDir:
        os.chdir(workDir)
        # else we must be there already
    os.chdir("./test")
   
    
    masterList = open("ReferenceDataset.txt")
    list = [l.strip('\n\r') for l in masterList.readlines()]
    return list


def bigLoop():
    referenceList = setup()
    suiteList = []
    for molecule in referenceList:
        sys.stderr.write(molecule + '\n')
        try:
            fileName = f"./{molecule}.suitegeom"
            inFile = open(fileName)
            suites = suitename.read(inFile)
            suites = suitename.compute(suites)
            for s in suites:
                s.source = molecule
                if s.issue:
                    s.situation = reasons[s.issue]
                # outSuite(s) # for debugging only
        except OSError as err:
            print(err)
        suiteList.append(suites)

    bigList = [item for sublist in suiteList for item in sublist]
    # suiteList is a list of lists; bigList is flattened
    bigList.sort(key=lambda suite: suite.comment + suite.situation) # swap + order to de-emphasize 7D
    for s in bigList:
        outSuite3(s)    # use 2 if de-emphasizing 7D


# modified especially for the sort order emphasizing 7D rejections
def outSuite3(suite):
    if not suite.valid: return
    
    note = ""
    reason = suite.comment + '  ' + suite.situation
    if suite.cluster.status == "wannabe":
        note = " wannabe"
    outIDs = ":".join(suite.pointID)
    output = (
        f"{suite.source}:{outIDs} {suite.bin.name} {suite.cluster.name} {float(suite.suiteness):5.3f}  {reason}{note}"
    )
    print(output)

# modified especially for the sort order
def outSuite2(suite):
    if not suite.valid: return
    
    note = ""
    reason = suite.situation + ' ' + suite.comment
    if suite.cluster.status == "wannabe":
        note = " wannabe"
    outIDs = ":".join(suite.pointID)
    output = (
        f"{suite.source}:{outIDs} {suite.bin.name} {suite.cluster.name} {float(suite.suiteness):5.3f}  {reason}{note}"
    )
    print(output)


# A fairly standard traditional output
def outSuite(suite):
    if not suite.valid: return

    reason = ""; note=""
    if suite.issue:
        reason = " " + reasons[suite.issue]
    elif suite.comment:
        reason = " " + suite.comment
    elif suite.situation:
        reason += " " + suite.situation
    if suite.cluster.status == "wannabe":
        note = " wannabe"
    outIDs = ":".join(suite.pointID)
    output = (
        f"{suite.source}:{outIDs} {suite.bin.name} {suite.cluster.name} {float(suite.suiteness):5.3f}{reason}{note}"
    )
    print(output)
bigLoop()

