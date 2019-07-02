from __future__ import absolute_import, division, print_function
from iotbx.bioinformatics import local_blast

#Note:
#I've used "X" to fill gaps and N-term positions to keep resnum consistant with
#blast output.  It seems to work fine. If you keep track of your sequences
#seperately, you probably can do it differently.   LWH 4/12/19

seq="XXXXNKLHVIDLHKRYGGHEVLKGVSLQARAGDVISIIGSSGSGKSTFLRCINFLEKPSEGAIIVNGQNINLVR\
DKDGQLKVADKNQLRLLRTRLTMVFQHFNLWSHMTVLENVMEAPIQVLGLSKHDARERALKYLAKVGIDERAQGKYPVH\
LSGGQQQRVSIARALAMEPDVLLFDEPTSALDPELVGEVLRIMQQLAEEGKTMVVVTHEMGFARHVSSHVIFLHQGKIE\
EEGDPEQVFGNPQSPRLQQFLKGSLKKLEH"

if __name__=="__main__":
  a=local_blast.pdbaa(seq=seq).run()
  xmldata="\n".join(a)
  hit="1B0U"
  assert hit in xmldata,"XML output not as expected. Pdbaa test failed."
  print("OK")
