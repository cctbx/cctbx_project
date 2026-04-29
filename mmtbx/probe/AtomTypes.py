##################################################################################
#                Copyright 2021  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License"],
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

##################################################################################
# This module exports a class that stores information about atoms that is needed by the
# Probe and Reduce portions of MolProbity to determine how to properly handle
# them.  This includes color information for Kinemage outputs but also information
# that may be available from CCTBX such as Van der Waals radii, donor/acceptor
# status, and whether the atom is metallic.
#
# Use the FindExtraAtomInfo() function to get an ExtraAtomInfo structure filled
# based on a specific atom along with a second tuple value that has a list of all
# information in the table for that atom.

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import, unicode_literals
import sys

##################################################################################
# Helper functions.

def Unpad(n):
  '''
  Gobble up all spaces from the end of the name.  Leave spaces at the beginning
  that come before non-space characters.
  '''
  while n[-1] == ' ':
    n = n[:-1]
  return n

# Is a carbon atom a Carbonyl from a standard amino acid?
# @todo This function is used once within Helpers.py to change the radius
# to match what is expectd based on experiments run by the Richardsons in
# September 2021.  If the CCTBX is changed to use the new values, then we
# will no longer need this function.
def IsSpecialAminoAcidCarbonyl(resName, atomName):
  """Given a residue and atom name, determine whether that atom is a C=O.
  This does not mark the ' C  ' atom that is always a Carbonyl; that is checked separately.
  :param resName: String containing the 1-3-character residue name in all caps, including leading space.
  :param atomName: String containing the 1-4-character atom name in all caps, including leading space.
  :returns True if the atom is a C=O in a standard residue, False if not.  Does not handle HET atoms.
  """
  if Unpad(atomName) == ' CG':
    return resName in ['ASP','ASN','ASX']
  if Unpad(atomName) == ' CD':
    return resName in ['GLU','GLN','GLX']
  return False

# Table of aromatic-acceptor atoms by residue and atom name.  The first entry in each list element is
# a list of residue names with trailing spaces trimmed.  The second is a list of atoms that qualify
# for all of the entries in the residue names.  In both cases, the strings are stripped of all
# spaces to the left and right.
_AromaticAcceptorTable = [
  # Note: Some atoms from these residues are listed in other sections.  The combination of
  # reside and atom name is not duplicated, but there are multiple entries for some residues --
  # this is not a set.

  # Do not mark the atoms in the Histidine ring or the atoms in the TRP 5-sided
  # ring as acceptors just because they are in the ring, but the HIS Nitrogens
  # remain acceptors because they might be unprotonated.  5-sided rings were getting too
  # many bonds from the sides.  Once we have a better aromatic ring hydrogen bond test,
  # we may put these back in.
  # @todo Remove the Hydrogens & Deuteriums below that are bonded to these atoms.
  [ ['HIS'], ['ND1','NE2'] ],
  # [ ['HIS'], ['CD2','CE1','CG'] ],
  # [ ['TRP'], ['CH2','CZ3','CZ2','CE3','CE2','NE1','CD2','CD1','CG'] ],
  [ ['TRP'], ['CH2','CZ3','CZ2','CE3','CE2','CD2'] ],

  [ ['ADE','A'], ['N1','N3','N7','C2','C4','C5','C6','C8','N9'] ],
  [ ['CYT','C'], ['N3','N1','C2','C4','C5','C6'] ],
  [ ['GUA','G'], ['N3','N7','N1','C2','C4','C5','C6','C8','N9'] ],
  [ ['THY','T'], ['N1','C2','N3','C4','C5','C6'] ],
  [ ['URA','U'], ['N1','C2','N3','C4','C5','C6'] ],

  [ ['DA'], ['N1','N3','N7','C2','C4','C5','C6','C8','N9'] ],
  [ ['DC'], ['N3','N1','C2','C4','C5','C6'] ],
  [ ['DG'], ['N3','N7','N1','C2','C4','C5','C6','C8','N9'] ],
  [ ['DT'], ['N1','C2','N3','C4','C5','C6'] ],

  [ ['HEM'], ['N A','N B','N C','N D'] ],

  # Here we treat the aromatic Pi-bonds as hydrogen bond acceptors.
  # Note: Some atoms from these residues are listed in other sections.  The combination of
  # reside and atom name is not duplicated, but there are multiple entries for some residues --
  # this is not a set.

  [ ['HEM'], ['C1A','C2A','C3A','C4A',
              'C1B','C2B','C3B','C4B',
              'C1C','C2C','C3C','C4C',
              'C1D','C2D','C3D','C4D'] ],
  [ ['PHE'], ['CZ','CE2','CE1','CD2','CD1','CG'] ],
  [ ['TYR'], ['CZ','CE2','CE1','CD2','CD1','CG'] ],
]

# Is a carbon or nitrogen or hydrogen atom part of an aromatic ring?
def IsAromaticAcceptor(resName, atomName):
  """Given a residue and atom name, determine whether that atom is and acceptor because it is
  part of an aromatic ring.
  :param resName: String containing the 1-3-character residue name in all caps, including leading space.
  :param atomName: String containing the 1-4-character atom name in all caps, including leading space.
  :returns True if the atom is aromatic in a standard residue, False if not.  Does not handle HET atoms.
  """

  for e in _AromaticAcceptorTable:
    if resName.strip() in e[0] and atomName.strip() in e[1]:
      return True
  return False


def Test(inFileName = None):

  #========================================================================
  # Check IsAromaticAcceptor() to ensure it gives results when expected and not when not.
  aromaticChecks = [
      ['PHE', 'CE2', True],
      ['  A',  'N1', True],
      ['ASN',   'O', False]
    ]
  for a in aromaticChecks:
    assert IsAromaticAcceptor(a[0],a[1]) == a[2], "AtomTypes.Test(): {} {} not marked as aromatic {}".format(a[0],a[1],a[2])

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB/CIF file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  # This will throw an assertion failure if there is a problem.
  Test(fileName)
  print('Success!')
