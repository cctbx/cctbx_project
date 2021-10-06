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
import re
import iotbx
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import mmtbx

import mmtbx_probe_ext as probe

##################################################################################
# Helper functions.

# Pad the name with spaces on the right to ensure that it is at least as long as
# requested.
def Pad(s, n=4):
  ret = s
  while len(s) < n:
    s += ' '
  return ret

# Gobble up all spaces from the end of the name after non-space characters
def Unpad(n):
  # Gobble up all spaces from the end of the name.
  while n[-1] == ' ':
    n = n[:-1]
  return n

# Is a carbon atom a Carbonyl from a standard amino acid?
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
_AromaticTable = [
  # Note: Some atoms from these residues are listed in other sections.  The combination of
  # reside and atom name is not duplicated, but there are multiple entries for some residues --
  # this is not a set.

  [ ['HIS'], ['ND1','NE2'] ],

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
  # [ ['HIS'], ['CD2','CE1','CG'] ],
  [ ['TRP'], ['CH2','CZ3','CZ2','CE3','CE2','NE1','CD2','CD1','CG'] ],

  # Here we add the hydrogens and deuteriums that can be part of a ring from probe:select.c
  [ ['PHE'], ['HD1','HD2','HE1','HE2','HZ','DD1','DD2','DE1','DE2','DZ'] ],
  [ ['HIS'], ['HD1','HD2','HE1','HE2','DD1','DD2','DE1','DE2'] ],
  [ ['TYR'], ['HD1','HD2','HE1','HE2','DD1','DD2','DE1','DE2'] ],
  [ ['TRP'], ['HD1','HE1','HE3','HZ2','HZ3','HH2','DD1','DE1','DE3','DZ2','DZ3','DH2'] ],
  [ ['U','URA','UTP','UDP','UMP','UR'], ['H3','HN3','H5','H6','D3','DN3','D5','D6'] ],
  [ ['T','THY','TTP','TDP','TMP','5MU','DT','TR'], ['H3','HN3','H6','D3','DN3','D6'] ],
  [ ['A','ADE','ATP','ADP','AMP','1MA','RIA','T6A','DA','AR'], ['H8','H2','D8','D2'] ],
  [ ['C','CYT','CTP','CDP','CMP','5MC','OMC','DC','CR'], ['H5','H6','D5','D6'] ],
  [ ['G','GUA','GTP','GDP','GMP','GSP','2MG','M2G','7MG','OMG','DG','GR'], ['H8','H1','HN1','D8','D1','DN1'] ],
  [ ['YG','1MG'], ['H8','D8'] ],
  [ ['PSU'], ['H6','D6','H1','HN1','D1','DN1','H3','HN3','D3','DN3'] ],
  [ ['I','DI'], ['H8','H2','H1','HN1','D8','D2','D1','DN1'] ]
]

# Is a carbon or nitrogen or hydrogen atom part of an aromatic ring?
def IsAromatic(resName, atomName):
  """Given a residue and atom name, determine whether that atom is part of an aromatic ring.
  :param resName: String containing the 1-3-character residue name in all caps, including leading space.
  :param atomName: String containing the 1-4-character atom name in all caps, including leading space.
  :returns True if the atom is aromatic in a standard residue, False if not.  Does not handle HET atoms.
  """

  for e in _AromaticTable:
    if resName.strip() in e[0] and atomName.strip() in e[1]:
      return True
  return False


##################################################################################

class AtomFlags(object):
  """Flags describing attributes that atoms can have.
  """
  EMPTY_FLAGS = 0               # No flags set
  IGNORE_ATOM = 1 << 0          # This atom should be ignored during processing, as if it did not exist
  DONOR_ATOM = 1 << 1           # Can be an electron donor
  ACCEPTOR_ATOM = 1 << 2        # Can be an electron acceptor
  HB_ONLY_DUMMY_ATOM = 1 << 3   # This is a dummy hydrogen added temporarily to a water when a donor is needed; it can Hbond but not clash.
  METALLIC_ATOM = 1 << 4        # The atom is metallic

##################################################################################

class AtomInfo(object):
  """Class that stores extra information about an atom that is looked up by the AtomTypes
  class methods.  The information is stored in properties.
  """

  def __init__(self, myValList = None):
    try:
      self._atomicNumber = myValList[0]
    except Exception:
      self._atomicNumber = 0
    try:
      self._name = myValList[1]
    except Exception:
      self._name = "?"
    try:
      self._fullName = myValList[2]
    except Exception:
      self._fullName = "unknown"
    try:
      self._vdwElectronCloudExplicit = myValList[3]
    except Exception:
      self._vdwElectronCloudExplicit = 0
    try:
      self._vdwNeutronExplicit = myValList[4]
    except Exception:
      self._vdwNeutronExplicit = 0
    try:
      self._vdwElectronCloudImplicit = myValList[5]
    except Exception:
      self._vdwElectronCloudImplicit = 0
    try:
      self._covalent = myValList[6]
    except Exception:
      self._covalent = 0
    try:
      self._kinemageColor = myValList[7]
    except Exception:
      self._kinemageColor = "grey"
    try:
      self._flags = myValList[8]
    except Exception:
      self._flags = AtomFlags.EMPTY_FLAGS

  # Getter and setter methods
  def get_atomicNumber(self): return self._atomicNumber
  def set_atomicNumber(self, val): self._atomicNumber = val
  def get_name(self): return self._name
  def set_name(self, val): self._name = val
  def get_fullName(self): return self._fullName
  def set_fullName(self, val): self._fullName = val
  def get_vdwElectronCloudExplicit(self): return self._vdwElectronCloudExplicit
  def set_vdwElectronCloudExplicit(self, val): self._vdwElectronCloudExplicit = val
  def get_vdwElectronCloudImplicit(self): return self._vdwElectronCloudImplicit
  def set_vdwElectronCloudImplicit(self, val): self._vdwElectronCloudImplicit = val
  def get_vdwNeutronExplicit(self): return self._vdwNeutronExplicit
  def set_vdwNeutronExplicit(self, val): self._vdwNeutronExplicit = val
  def get_covalent(self): return self._covalent
  def set_covalent(self, val): self._covalent = val
  def get_kinemageColor(self): return self._kinemageColor
  def set_kinemageColor(self, val): self._kinemageColor = val
  def get_flags(self): return self._flags
  def set_flags(self, val): self._flags = val

  # Properties
  atomicNumber = property(get_atomicNumber, set_atomicNumber)
  name = property(get_name, set_name)
  fullName = property(get_fullName, set_fullName)
  vdwElectronCloudExplicit = property(get_vdwElectronCloudExplicit, set_vdwElectronCloudExplicit)
  vdwElectronCloudImplicit = property(get_vdwElectronCloudImplicit, set_vdwElectronCloudImplicit)
  vdwNeutronExplicit = property(get_vdwNeutronExplicit, set_vdwNeutronExplicit)
  covalent = property(get_covalent, set_covalent)
  kinemageColor = property(get_kinemageColor, set_kinemageColor)
  flags = property(get_flags, set_flags)

class AtomTypes(object):
  """Class that looks up extra information for atoms that is required by the MolProbity Probe and
  Reduce modules.
  """

  def __init__(self, useNeutronDistances = False, useImplicitHydrogenDistances = False):
    """Constructor.
    :param useNeutronDistances: Use neutron distances and radii for scoring.
    The default is to use electron-cloud distances.  This is used both for the
    separation between a Hydgrogen and its bound partner and for the radius of the
    Hydrogen and it must be set consistently across the entire code base.  When this is
    True, it supercedes useImplicitHydrogenDistances.
    :param useImplicitHydrogenDistances: Default is to use distances consistent with
    explicitly-listed Hydrgoens, but setting this to True implicit-Hydrogen distances instead.
    This must be set consistently with the hydrogens in the model.
    """

    ##################################################################################
    # Store state based on options.
    self._useNeutronDistances = useNeutronDistances
    self._useImplicitHydrogenDistances = useImplicitHydrogenDistances

    ##################################################################################
    # Table of information about each type of atom.  The elements in each list are as
    # follows:
    #   Atomic number
    #   Name of the type, used to look up the atom type
    #   Full name of the type, useful when printing
    #   VDW radius for explicit hydrogen bonds at the electron cloud distance
    #   VDW radius for explicit hydrogen bonds at the neutron distance
    #   VDW radius for implicit hydrogen bonds at the electron cloud distance
    #   Covalent bond radius
    #   Name of the color to use in Mage/Kinemage to display the atom
    #   Flags describing the behavior of the atom, as follows:
    #     IGNORE_ATOM : This atom should be ignored
    #     DONOR_ATOM : This atom can be a hydrogen bond donor
    #     ACCEPTOR_ATOM : This atom can be a hydrogen bond acceptor
    #     HB_ONLY_DUMMY_ATOM : This is a dummy hydrogen added temporarily
    #         to a water when a donor is needed; it can Hbond but not clash.
    #     METALLIC_ATOM : This atom is metallic
    #
    # This table is based on the following:
    #   For non-metals, explicit VDW radii from
    #   Gavezzotti, J. Am. Chem. Soc. (1983) 105, 5220-5225.
    #   or, if unavailable,
    #   Bondi, J. Phys. Chem. (1964), V68, N3, 441-451.
    #   Covalent and ionic radii from
    #   Advanced Inorganic Chemistry, Cotton & Wilkinson, 1962, p93.

    self._AtomTable = [
      [ 0, "?",  "unknown",            1.05, 1.05, 0.00, 0.00, "magenta", AtomFlags.EMPTY_FLAGS],
      [ 0, "ignore", "ignore",         0.00, 0.00, 0.00, 0.00, "magenta", AtomFlags.IGNORE_ATOM],
      [ 1, "H",  "hydrogen",           1.22, 1.17, 0.00, 0.30, "grey",   AtomFlags.EMPTY_FLAGS],
      [ 1, "Har","hydrogen(aromatic)", 1.05, 1.00, 0.00, 0.30, "grey",   AtomFlags.EMPTY_FLAGS],
      [ 1, "Hpol","hydrogen(polar)",   1.05, 1.00, 0.00, 0.30, "grey",   AtomFlags.DONOR_ATOM],
      [ 1, "Ha+p","hydrogen(aromatic&polar)", 1.05, 1.00, 0.00, 0.30, "grey",   AtomFlags.DONOR_ATOM],
      [ 1, "HOd","hydrogen(only dummy)", 1.05, 1.00, 0.00, 0.30, "grey",   AtomFlags.DONOR_ATOM|AtomFlags.HB_ONLY_DUMMY_ATOM],
      [ 6, "C",  "carbon",             1.70, 1.70, 1.90, 0.77, "white",  AtomFlags.EMPTY_FLAGS],
      [ 6, "Car","carbon(aromatic)",   1.75, 1.75, 1.90, 0.77, "white",  AtomFlags.ACCEPTOR_ATOM],
      [ 6, "C=O","carbon(carbonyl)",   1.65, 1.65, 1.80, 0.77, "white",  AtomFlags.EMPTY_FLAGS],
      [ 7, "N",  "nitrogen",           1.55, 1.55, 1.70, 0.70, "sky",    AtomFlags.EMPTY_FLAGS],
      [ 7, "Nacc","nitrogen(acceptor)",1.55, 1.55, 1.70, 0.70, "sky",    AtomFlags.ACCEPTOR_ATOM],
      [ 8, "O",  "oxygen",             1.40, 1.40, 1.50, 0.66, "red",    AtomFlags.ACCEPTOR_ATOM],
      [15, "P",  "phosphorus",         1.80, 1.80, 1.80, 1.10, "pink",   AtomFlags.EMPTY_FLAGS],
      [16, "S",  "sulfur",             1.80, 1.80, 1.90, 1.04, "yellow", AtomFlags.ACCEPTOR_ATOM],
      [33, "As", "arsenic",            2.00, 2.00, 2.10, 1.21, "grey",   AtomFlags.EMPTY_FLAGS],
      [34, "Se", "selenium",           1.90, 1.90, 2.00, 1.17, "green",  AtomFlags.EMPTY_FLAGS],
      [ 9, "F",  "fluorine",           1.30, 1.30, 1.30, 0.58, "green",  AtomFlags.ACCEPTOR_ATOM],
      [17, "Cl", "chlorine",           1.77, 1.77, 1.77, 0.99, "green",  AtomFlags.ACCEPTOR_ATOM],
      [35, "Br", "bromine",            1.95, 1.95, 1.95, 1.14, "brown",  AtomFlags.ACCEPTOR_ATOM],
      [53, "I",  "iodine",             2.10, 2.10, 2.10, 1.33, "brown",  AtomFlags.ACCEPTOR_ATOM],

      # for most common metals we use Pauling's ionic radii
      # "covalent radii" does not really mean anything, but the interaction distance for
      # all interactions uses the same, ionic, radius.  This code differs from the C++
      # and C tables from Probe and Reduce.  It was modified because of further study by the
      # Richardson group in 2021.
      [ 3, "Li", "lithium",            0.60, 0.60, 0.60, 0.60, "grey", AtomFlags.METALLIC_ATOM],
      [11, "Na", "sodium",             0.95, 0.95, 0.95, 0.95, "grey", AtomFlags.METALLIC_ATOM],
      [13, "Al", "aluminum",           0.50, 0.50, 0.50, 0.50, "grey", AtomFlags.METALLIC_ATOM],
      [19, "K",  "potassium",          1.33, 1.33, 1.33, 1.33, "grey", AtomFlags.METALLIC_ATOM],
      [12, "Mg", "magnesium",          0.65, 0.65, 0.65, 0.65, "grey", AtomFlags.METALLIC_ATOM],
      [20, "Ca", "calcium",            0.99, 0.99, 0.99, 0.99, "grey", AtomFlags.METALLIC_ATOM],
      [25, "Mn", "manganese",          0.80, 0.80, 0.80, 0.80, "grey", AtomFlags.METALLIC_ATOM],
      [26, "Fe", "iron",               0.74, 0.74, 0.74, 0.74, "grey", AtomFlags.METALLIC_ATOM],
      [27, "Co", "cobalt",             0.70, 0.70, 0.70, 0.70, "blue", AtomFlags.METALLIC_ATOM],
      [28, "Ni", "nickel",             0.66, 0.66, 0.66, 0.66, "grey", AtomFlags.METALLIC_ATOM],
      [29, "Cu", "copper",             0.72, 0.72, 0.72, 0.72,"orange",AtomFlags.METALLIC_ATOM],
      [30, "Zn", "zinc",               0.71, 0.71, 0.71, 0.71, "grey", AtomFlags.METALLIC_ATOM],
      [37, "Rb", "rubidium",           1.48, 1.48, 1.48, 1.48, "grey", AtomFlags.METALLIC_ATOM],
      [38, "Sr", "strontium",          1.10, 1.10, 1.10, 1.10, "grey", AtomFlags.METALLIC_ATOM],
      [42, "Mo", "molybdenum",         0.93, 0.93, 0.93, 0.93, "grey", AtomFlags.METALLIC_ATOM],
      [47, "Ag", "silver",             1.26, 1.26, 1.26, 1.26, "white",AtomFlags.METALLIC_ATOM],
      [48, "Cd", "cadmium",            0.91, 0.91, 0.91, 0.91, "grey", AtomFlags.METALLIC_ATOM],
      [49, "In", "indium",             0.81, 0.81, 0.81, 0.81, "grey", AtomFlags.METALLIC_ATOM],
      [55, "Cs", "cesium",             1.69, 1.69, 1.69, 1.69, "grey", AtomFlags.METALLIC_ATOM],
      [56, "Ba", "barium",             1.29, 1.29, 1.29, 1.29, "grey", AtomFlags.METALLIC_ATOM],
      [79, "Au", "gold",               1.10, 1.10, 1.10, 1.10, "gold", AtomFlags.METALLIC_ATOM],
      [80, "Hg", "mercury",            1.00, 1.00, 1.00, 1.00, "grey", AtomFlags.METALLIC_ATOM],
      [81, "Tl", "thallium",           1.44, 1.44, 1.44, 1.44, "grey", AtomFlags.METALLIC_ATOM],
      [82, "Pb", "lead",               0.84, 0.84, 0.84, 0.84, "grey", AtomFlags.METALLIC_ATOM],

      # for other metals we use Shannon's ionic radii
      # Acta Crystallogr. (1975) A32, pg751.
      [23, "V",  "vanadium",           0.79, 0.79, 0.79, 0.79, "grey", AtomFlags.METALLIC_ATOM],
      [24, "Cr", "chromium",           0.73, 0.73, 0.73, 0.73, "grey", AtomFlags.METALLIC_ATOM],
      [52, "Te", "tellurium",          0.97, 0.97, 0.97, 0.97, "grey", AtomFlags.METALLIC_ATOM],
      [62, "Sm", "samarium",           1.08, 1.08, 1.08, 1.08, "grey", AtomFlags.METALLIC_ATOM],
      [64, "Gd", "gadolinium",         1.05, 1.05, 1.05, 1.05, "grey", AtomFlags.METALLIC_ATOM],
      [70, "Yb", "ytterbium",          1.14, 1.14, 1.14, 1.14, "grey", AtomFlags.METALLIC_ATOM],
      [74, "W",  "tungsten",           0.66, 0.66, 0.66, 0.66, "grey", AtomFlags.METALLIC_ATOM],
      [78, "Pt", "platinum",           0.63, 0.63, 0.63, 0.63, "grey", AtomFlags.METALLIC_ATOM],
      [92, "U",  "uranium",            1.03, 1.03, 1.03, 1.03, "grey", AtomFlags.METALLIC_ATOM],

      # Cotton & Wilkinson and also-
      # L.E. Sutton (ed.) in Table of interatomic distances and configuration in molecules
      # and ions, Supplement 1956-1959, Special publication No. 18, Chemical Society,
      # London, UK, 1965 (as listed in web-elements by Mark Winter)
      # http://www.shef.ac.uk/chemistry/web-elements
      [ 2, "He",  "helium",            1.60, 1.60, 1.60, 0.00, "sky",          AtomFlags.EMPTY_FLAGS],
      [ 4, "Be",  "beryllium",         0.31, 0.31, 0.31, 0.90, "grey", AtomFlags.METALLIC_ATOM],
      [ 5, "B",   "boron",             0.20, 0.20, 0.20, 0.86, "grey",         AtomFlags.EMPTY_FLAGS],
      [10, "Ne",  "neon",              1.60, 1.60, 1.60, 0.00, "pink",         AtomFlags.EMPTY_FLAGS],
      [14, "Si",  "silicon",           2.10, 2.10, 2.10, 1.17, "grey", AtomFlags.METALLIC_ATOM],
      [18, "Ar",  "argon",             1.89, 1.89, 1.89, 0.00, "orange",       AtomFlags.EMPTY_FLAGS],
      [21, "Sc",  "scandium",          0.68, 0.68, 0.68, 0.44, "grey", AtomFlags.METALLIC_ATOM],
      [22, "Ti",  "titanium",          0.75, 0.75, 0.75, 1.49, "grey", AtomFlags.METALLIC_ATOM],
      [31, "Ga",  "gallium",           0.53, 0.53, 0.53, 1.27, "grey", AtomFlags.METALLIC_ATOM],
      [32, "Ge",  "germanium",         0.60, 0.60, 0.60, 1.34, "grey", AtomFlags.METALLIC_ATOM],
      [36, "Kr",  "krypton",           2.01, 2.01, 2.01, 0.00, "greentint",    AtomFlags.EMPTY_FLAGS],
      [39, "Y",   "yttrium",           0.90, 0.90, 0.90, 1.64, "grey", AtomFlags.METALLIC_ATOM],
      [40, "Zr",  "zirconium",         0.77, 0.77, 0.77, 1.51, "grey", AtomFlags.METALLIC_ATOM],
      [50, "Sn",  "tin",               0.71, 0.71, 0.71, 1.45, "grey", AtomFlags.METALLIC_ATOM],
      [51, "Sb",  "antimony",          2.20, 2.20, 2.20, 1.41, "grey", AtomFlags.METALLIC_ATOM],
      [54, "Xe",  "xenon",             2.18, 2.18, 2.18, 0.00, "magenta",      AtomFlags.EMPTY_FLAGS],
      [57, "La",  "lanthanum",         1.03, 1.03, 1.03, 1.77, "grey", AtomFlags.METALLIC_ATOM],
      [58, "Ce",  "cerium",            0.87, 0.87, 0.87, 1.61, "grey", AtomFlags.METALLIC_ATOM],
      [87, "Fr",  "francium",          1.94, 1.94, 1.94, 2.68, "grey", AtomFlags.METALLIC_ATOM],
      [88, "Ra",  "radium",            1.62, 1.62, 1.62, 2.36, "grey", AtomFlags.METALLIC_ATOM],
      [90, "Th",  "thorium",           1.08, 1.08, 1.08, 1.82, "grey", AtomFlags.METALLIC_ATOM],

      # finally, we have a set of elements where the radii are unknown
      # so we use estimates and extrapolations based on web-elements data
      [41, "Nb",  "niobium",           0.86, 0.86, 0.86, 1.40, "grey", AtomFlags.METALLIC_ATOM],
      [43, "Tc",  "technetium",        0.71, 0.71, 0.71, 1.25, "grey", AtomFlags.METALLIC_ATOM],
      [44, "Ru",  "ruthenium",         0.82, 0.82, 0.82, 1.36, "grey", AtomFlags.METALLIC_ATOM],
      [45, "Rh",  "rhodium",           0.76, 0.76, 1.76, 1.30, "grey", AtomFlags.METALLIC_ATOM],
      [46, "Pd",  "palladium",         1.05, 1.05, 1.05, 1.59, "grey", AtomFlags.METALLIC_ATOM],
      [59, "Pr",  "praseodymium",      1.11, 1.11, 1.11, 1.65, "grey", AtomFlags.METALLIC_ATOM],
      [60, "Nd",  "neodymium",         1.10, 1.10, 1.10, 1.64, "grey", AtomFlags.METALLIC_ATOM],
      [61, "Pm",  "promethium",        1.15, 1.15, 1.15, 1.89, "grey", AtomFlags.METALLIC_ATOM],
      [63, "Eu",  "europium",          1.31, 1.31, 1.31, 1.85, "grey", AtomFlags.METALLIC_ATOM],
      [65, "Tb",  "terbium",           1.05, 1.05, 1.05, 1.59, "grey", AtomFlags.METALLIC_ATOM],
      [66, "Dy",  "dysprosium",        1.05, 1.05, 1.05, 1.59, "grey", AtomFlags.METALLIC_ATOM],
      [67, "Ho",  "holmium",           1.04, 1.04, 1.04, 1.58, "grey", AtomFlags.METALLIC_ATOM],
      [68, "Er",  "erbium",            1.03, 1.03, 1.03, 1.57, "grey", AtomFlags.METALLIC_ATOM],
      [69, "Tm",  "thulium",           1.02, 1.02, 1.02, 1.56, "grey", AtomFlags.METALLIC_ATOM],
      [71, "Lu",  "lutetium",          1.02, 1.02, 1.02, 1.56, "grey", AtomFlags.METALLIC_ATOM],
      [72, "Hf",  "hafnium",           0.85, 0.85, 0.85, 1.46, "grey", AtomFlags.METALLIC_ATOM],
      [73, "Ta",  "tantalum",          0.86, 0.86, 0.86, 1.40, "grey", AtomFlags.METALLIC_ATOM],
      [75, "Re",  "rhenium",           0.77, 0.77, 0.77, 1.31, "grey", AtomFlags.METALLIC_ATOM],
      [76, "Os",  "osmium",            0.78, 0.78, 0.78, 1.32, "grey", AtomFlags.METALLIC_ATOM],
      [77, "Ir",  "iridium",           0.80, 0.80, 0.80, 1.34, "grey", AtomFlags.METALLIC_ATOM],
      [83, "Bi",  "bismuth",           1.17, 1.17, 1.17, 1.71, "grey", AtomFlags.METALLIC_ATOM],
      [84, "Po",  "polonium",          0.99, 0.99, 0.99, 1.53, "grey", AtomFlags.METALLIC_ATOM],
      [85, "At",  "astatine",          0.91, 0.91, 0.91, 1.45, "grey", AtomFlags.METALLIC_ATOM],
      [86, "Rn",  "radon",             2.50, 2.50, 2.50, 0.00, "pinktint",     AtomFlags.EMPTY_FLAGS],
      [89, "Ac",  "actinium",          1.30, 1.30, 1.30, 2.00, "grey", AtomFlags.METALLIC_ATOM],
      [91, "Pa",  "protoactinium",     1.10, 1.10, 1.10, 1.85, "grey", AtomFlags.METALLIC_ATOM],
      [93, "Np",  "neptunium",         1.00, 1.00, 1.00, 1.72, "grey", AtomFlags.METALLIC_ATOM],
      [94, "Pu",  "plutonium",         1.00, 1.00, 1.00, 1.67, "grey", AtomFlags.METALLIC_ATOM],
      [95, "Am",  "americium",         1.00, 1.00, 1.00, 1.63, "grey", AtomFlags.METALLIC_ATOM],
      [96, "Cm",  "curium",            1.00, 1.00, 1.00, 1.60, "grey", AtomFlags.METALLIC_ATOM],
      [97, "Bk",  "berkelium",         1.00, 1.00, 1.00, 1.58, "grey", AtomFlags.METALLIC_ATOM],
      [98, "Cf",  "californium",       1.00, 1.00, 1.00, 1.57, "grey", AtomFlags.METALLIC_ATOM],
      [99, "Es",  "einsteinium",       1.00, 1.00, 1.00, 1.56, "grey", AtomFlags.METALLIC_ATOM],
      [100,"Fm",  "fermium",           1.00, 1.00, 1.00, 1.55, "grey", AtomFlags.METALLIC_ATOM],
      [101,"Md",  "mendelevium",       1.00, 1.00, 1.00, 1.55, "grey", AtomFlags.METALLIC_ATOM],
      [102,"No",  "nobelium",          1.00, 1.00, 1.00, 1.55, "grey", AtomFlags.METALLIC_ATOM]
    ]

    ##################################################################################
    # Construct a dictionary that maps from the name of the type in the _AtomTable
    # to its full entry in the table to make it fast to look up an atom by its type
    # name.
    self._Index = {}
    for e in self._AtomTable:
      self._Index[e[1]] = e

    ##################################################################################
    # Make a string that has all of the special first characters from an atom name that
    # would cause it to parse the remainder of the name as if it were a full name (this
    # uses a table that is a subset of the full name-parsing table).
    self._specialAtomFirstChars = '*"\'`_+- 0123456789'

    ##################################################################################
    # Make a dictionary for upper-case second letters in multi-letter atom names where
    # the first letter is H but the atom is not a Hydrogen.
    # It stores the list of residues for which this name is valid.
    # For example, "Hg" would have an entry named 'G' which lists the residues that can
    # have an atom named "Hg" in it.
    # Making an empty list here will have the same effect as not having an entry for
    # a given letter.  All names in the list must be fully upper case.  Spaces are
    # significant in the residue names.
    # Atoms whose names start with H but are not on one of the lists will be converted
    # to Hydrogen.
    self._legalResiduesForHElements = {
      'E' : [],
      'F' : ['PHF', 'HF3', 'HF5'],
      'G' : [' HG', 'HG2', 'HGB', 'HGC', 'HGI', 'MAC', 'MBO', 'MMC',
             'PHG', 'PMB', 'AAS', 'AMS', 'BE7', 'CMH', 'EMC', 'EMT'],
      'O' : [' HO', 'HO3'],
      'S' : []
    }

    ##################################################################################
    # Upper-cased legal names for atoms whose name begins with H.  Used to generate a
    # warning if we have a different name than these and it is not a special atom name.
    # HH and HD are there to get rid of warnings for PBD v3 names and H5'' is
    # to handle RNA/DNA backbone atoms.
    self._legalStandardHAtomNames = ['H', 'HE', 'HF', 'HG', 'HO', 'HS',
      "H5''", 'HH', 'HD']

    ##################################################################################
    # Table of allowed names for atoms that begin with one of the special characters.
    # This table operates on the name after the special character has been removed.
    # Each entry has a regular expression to match, the resulting name, and a Boolean
    # telling whether to warn about this translation.
    self._specialNameTable = [
      [ r'A.1',    'O',  True ],
      [ r'A.2',    'N',  True ],
      [ r'B.*',    'B',  False ],
      [ r'C.*',    'C',  False ],
      [ r'D.*',    'H',  True ],  # These are counted as H internally, but output at D
      [ r'F.*',    'F',  False ],
      # H atoms are handled separately
      [ r'I.*',    'I',  False ],
      [ r'K.*',    'K',  False ],
      [ r'N.*',    'N',  False ],
      [ r'O.*',    'O',  False ],
      [ r'P.*',    'P',  False ],
      [ r'S.*',    'S',  False ],
      [ r'U.*',    'U',  False ],
      [ r'V.*',    'V',  False ],
      [ r'W.*',    'W',  False ],
      [ r'Y.*',    'Y',  False ]
    ]

    ##################################################################################
    # Table of allowed names for atoms that do not begin with one of the special characters.
    # Each entry has a regular expression to match, the resulting name, and a Boolean
    # telling whether to warn about this translation.
    # The dot character '.' matches any single character and the * character means 0 or more
    # instances of the previous character, so '.*' at the end matches an arbitrary ending.
    # The expressions should be upper-case.  The resulting names should match the case of
    # the _Index dictionary.
    # Earlier entries are checked before later entries, so for example the C.* at the end of the
    # C's will match all atoms that don't match other than the C entries before it.
    self._nameTable = [
      [ r'AC.*',    'C',  True ],
      [ r'AG.*',    'Ag', False ],
      [ r'AH.*',    'H',  True ],
      [ r'AL.*',    'Al', False ],
      [ r'AM.*',    'Am', False ],
      [ r'AN.*',    'N',  True ],
      [ r'AO.*',    'O',  True ],
      [ r'AP.*',    'P',  True ],
      [ r'AR.*',    'Ar', False ],
      [ r'AS.*',    'As', False ],
      [ r'AT.*',    'At', False ],
      [ r'AU.*',    'Au', False ],

      [ r'BA.*',    'Ba', False ],
      [ r'BE.*',    'Be', False ],
      [ r'BI.*',    'Bi', False ],
      [ r'BK.*',    'Bk', False ],
      [ r'BR.*',    'Br', False ],

      [ r'CA.*',    'Ca', False ],
      [ r'CC.*',    'C',  True ],
      [ r'CD.*',    'Cd', False ],
      [ r'CE.*',    'Ce', False ],
      [ r'CF.*',    'Cf', False ],
      [ r'CH.*',    'H',  True ],
      [ r'CL.*',    'Cl', False ],
      [ r'CM.*',    'Cm', False ],
      [ r'CN.*',    'N',  True ],
      [ r'CO.*',    'Co', False ],
      [ r'CP.*',    'C',  True ],
      [ r'CR.*',    'Cr', False ],
      [ r'CS.*',    'Cs', False ],
      [ r'CU.*',    'Cu', False ],
      [ r'C.*',     'C',  True ],  # All other atoms starting with C are called Carbon with a warning

      [ r'DY.*',    'Dy', False ],
      [ r'DC.*',    'C',  True ],
      [ r'DH.*',    'H',  True ],
      [ r'DN.*',    'N',  True ],
      [ r'DO.*',    'O',  True ],
      [ r'DP.*',    'P',  True ],
      [ r'D.*',     'H',  True ],   # All other atoms starting with D are called Hydrogen

      [ r'ER.*',    'Er', False ],
      [ r'ES.*',    'Es', False ],
      [ r'EU.*',    'Eu', False ],
      [ r'EC.*',    'C',  True ],
      [ r'EH.*',    'H',  True ],
      [ r'EN.*',    'N',  True ],
      [ r'EO.*',    'O',  True ],
      [ r'EP.*',    'P',  True ],

      [ r'FE.*',    'Fe', False ],
      [ r'FM.*',    'Fm', False ],
      [ r'FR.*',    'Fr', False ],
      [ r'FC.*',    'C',  True ],
      [ r'FH.*',    'H',  True ],
      [ r'FN.*',    'N',  True ],
      [ r'FO.*',    'O',  True ],
      [ r'FP.*',    'P',  True ],

      [ r'GA.*',    'Ga', False ],
      [ r'GD.*',    'Gd', False ],
      [ r'GE.*',    'Ge', False ],
      [ r'GC.*',    'C',  True ],
      [ r'GH.*',    'H',  True ],
      [ r'GN.*',    'N',  True ],
      [ r'GO.*',    'O',  True ],
      [ r'GP.*',    'P',  True ],

      # H atoms are handled separately

      [ r'IN.*',    'In', False ],
      [ r'IR.*',    'Ir', False ],

      [ r'KR.*',    'Kr', False ],

      [ r'LA.*',    'La', False ],
      [ r'LI.*',    'Li', False ],
      [ r'LU.*',    'Lu', False ],

      [ r'MD.*',    'Md', False ],
      [ r'MG.*',    'Mg', False ],
      [ r'MN.*',    'Mn', False ],
      [ r'MO.*',    'Mo', False ],

      [ r'NA.*',    'Na', True ],
      [ r'NB.*',    'Nb', True ],
      [ r'NC.*',    'C',  True ],
      [ r'ND.*',    'Nd', True ],
      [ r'NE.*',    'Ne', True ],
      [ r'NH.*',    'H',  True ],
      [ r'NI.*',    'Ni', False ],
      [ r'NN.*',    'N',  True ],
      [ r'NO.*',    'O',  True ],   # Non standard
      [ r'NP.*',    'P',  True ],   # Non standard
      [ r'NS.*',    'S',  True ],
      [ r'N.*',     'N',  True ],   # All other atoms starting with N are called Nitrogen

      [ r'OS.*',    'Os', False ],
      [ r'O.*',     'O',  True ],   # All other atoms starting with O are called Oxygen

      [ r'PA.*',    'Pa', True ],
      [ r'PB.*',    'Pb', True ],
      [ r'PD.*',    'Pd', True ],
      [ r'PM.*',    'Pm', False ],
      [ r'PO.*',    'Po', False ],
      [ r'PR.*',    'Pr', False ],
      [ r'PT.*',    'Pt', False ],
      [ r'PU.*',    'Pu', False ],
      [ r'P.*',     'P',  True ],   # All other atoms starting with P

      [ r'RA.*',    'Ra', False ],
      [ r'RB.*',    'Rb', False ],
      [ r'RE.*',    'Re', False ],
      [ r'RH.*',    'Rh', False ],
      [ r'RN.*',    'Rn', False ],
      [ r'RU.*',    'Ru', False ],

      [ r'SB.*',    'Sb', True ],
      [ r'SC.*',    'Sc', False ],
      [ r'SE.*',    'Se', True ],
      [ r'SI.*',    'Si', False ],
      [ r'SM.*',    'Sm', False ],
      [ r'SN.*',    'Sn', False ],
      [ r'SR.*',    'Sr', False ],
      [ r'S.*',     'S',  True ],   # All other atoms starting with S

      [ r'TA.*',    'Ta', False ],
      [ r'TB.*',    'Tb', False ],
      [ r'TC.*',    'Tc', False ],
      [ r'TE.*',    'Te', False ],
      [ r'TH.*',    'Th', False ],
      [ r'TI.*',    'Ti', False ],
      [ r'TL.*',    'Tl', False ],
      [ r'TM.*',    'Tm', False ],

      [ r'XE.*',    'Xe', False ],

      [ r'YB.*',    'Yb', False ],

      [ r'ZN.*',    'Zn', False ],
      [ r'ZR.*',    'Zr', False ]
    ]

    ##################################################################################
    # Table of last-chance names for atoms that were not found in one of the above
    # tables.
    # This table operates on the name after the first character has been removed.
    # Each entry has a regular expression to match, the resulting name, and a Boolean
    # telling whether to warn about this translation -- these always warn.
    self._lastChanceNameTable = [
      [ r'H.*',    'H',  True ],
      [ r'D.*',    'H',  True ],
      [ r'C.*',    'C',  True ],
      [ r'N.*',    'N',  True ],
      [ r'O.*',    'O',  True ],
      [ r'P.*',    'P',  True ],
      [ r'S.*',    'S',  True ],

      [ r'I.*',    'I',  True ],
      [ r'K.*',    'K',  True ],
      [ r'V.*',    'V',  True ],
      [ r'W.*',    'W',  True ],
      [ r'U.*',    'U',  True ],

      [ r'AG.*',   'Ag', True ],
      [ r'AL.*',   'Al', True ],
      [ r'AS.*',   'As', True ],
      [ r'AU.*',   'Au', True ],

      [ r'FE.*',   'Fe', True ],

      [ r'GD.*',   'Gd', True ],

      [ r'LI.*',   'Li', True ],

      [ r'MG.*',   'Mg', True ],
      [ r'MN.*',   'Mn', True ],
      [ r'MO.*',   'Mo', True ],

      [ r'ZN.*',   'Zn', True ]
    ]

  ##################################################################################
  # Given an iotbx.pdb.atom, look up its mmtbx_probe_ext.ExtraAtomInfo in the atom table.
  # This includes checking for special cases using the helper functions above.

  def FindAtomInfo(self, atom):
    """Given an iotbx.pdb.atom, look up its information in the atom table.
    :param atom: iotbx.pdb.atom entry to look up.
    :returns a pair (filled in AtomInfo, string warning) on success, raises ValueError on failure.
    The warning will be empty if there was no problem with the lookup and will be a
    printable string explaining the problem if there was one.
    """

    # The element name we're going to use to look up in the table.
    elementName = None

    # Should we emit a warning about this name translation?
    emitWarning = False

    # Find the name of the atom and the residue it is in.  The atom's parent is an
    # atom group, which holds its residue name.  The atom name is padded on the right
    # with spaces so that it is at least four characters long and made upper-case.
    atomPadded = Pad(atom.name.upper())
    atomName = atomPadded # name that will be adjusted as we go
    resName = atom.parent().resname.upper()

    # Special-case check for a mis-placed Selenium atom
    if atomName[0:2] == ' SE':
      elementName = 'Se'
      emitWarning = True

    # Based on fixupAmbigAtomName() from Reduce ElementInfo.cpp
    # Fix up ambiguous atom names before trying to figure out what they are.  Emit a
    # warning for each of them if we do rename them.
    if resName in ['ASN','GLN'] and atomName in [' AD1',' AD2',' AE1',' AE2']:
      emitWarning = True
      if atomName[3] == '1':
        atomName[1] = 'O'
      else:
        # All of the other entries in the list have '2', so we don't need a separate check
        atomName[1] = 'N'

    ###
    # If we didn't hit a special case above, check using the appropriate table
    # based on the first character of the element name.
    if elementName is None:
      # See if the first character of the atom name is special.  If so, shift it off
      # the name and replace the standard parsing table with the subset table for special
      # atom names.
      nameTable = self._nameTable
      specialName = False
      if atomName[0] in self._specialAtomFirstChars:
        atomName = atomName[1:]
        nameTable = self._specialNameTable
        specialName = True

      # For Hydrogens, see if we need to truncate the name by removing its second character
      # because it is an He, Hf, Hg, Ho, or Hs in a residue that does not allow such names.
      # If it is disallowed, replace it with a simple 'H'.
      if atomName[0] == 'H':
        # If we are not a truncated special name, then we do more scrutiny on the name and
        # emit a warning if it is not any of the recognized names.
        if (not specialName) and (not atomName in self._legalStandardHAtomNames):
          emitWarning = True

        # Use the default value of 'H' unless the second letter is in a list of residues
        # that allow that second letter.
        elementName = 'H'
        try:
          if resName in self._legalResiduesForHElements[atomName[1]]:
            # This residue is in the list of those that are valid for this name,
            # so we make the element name H followed by the same lower-case letter.
            elementName = 'H' + atomName[1].lower()
        except KeyError:
          # We did not find an entry for this character, so we leave things alone
          pass

      # Look up the atom in all entries of the table to see if any of their regular
      # expressions match.  If so, set the elementName and warning emission based on
      # the table entry.
      for n in nameTable:
        e = re.compile(n[0])
        if e.match(atomName) is not None:
          elementName = n[1]
          emitWarning |= n[2]
          break;

    ###
    # If we did not find an elementName yet, we always emit a warning,
    # then default to Carbon, and then try another pass on names
    # always skipping the first character of the name.
    if elementName is None:
      elementName = 'C'
      emitWarning = True

      atomName = atomPadded[1:]
      for n in self._lastChanceNameTable:
        e = re.compile(n[0])
        if e.match(atomName) is not None:
          elementName = n[1]
          emitWarning |= n[2]
          break;

    ###
    # Some Carbon atoms are too general when coming from the tables, so we need to adjust
    # them to match the correct types.
    # Based on reduce.cpp fixupHeavyAtomElementType().

    # This maps from amino acid names to single letter codes, but we just care if it is in the list.
    aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
    if elementName == 'C':
      if IsAromatic(resName, atomPadded):
        elementName = 'Car'
      elif atomPadded == ' C  ' and atom.parent().resname in aa_resnames:
        elementName = 'C=O'
      elif IsSpecialAminoAcidCarbonyl(resName, atomPadded):
        elementName = 'C=O'
      # @todo Aromatic and carbonyl carbons in non-standard HET residues need to be identified somehow
      # (C++ code had this comment in it without code to handle the situation).

    ###
    # Add Acceptor flag for some Nitrogen atoms.
    # Based on reduce.cpp fixupHeavyAtomElementType().
    if elementName == 'N':
      # All acceptors in the C++ table were also marked as aromatic and vice-versa,
      # so we use the test for both to locate them.
      if IsAromatic(resName, atomPadded):
        elementName = 'Nacc'
      #elif @todo Het N atoms with 1 or 2 bonded (connected) neighbors:
      #  elementName = 'Nacc'
      #elif @todo This is a fragment and may be an N-terminal (reduce.cpp line 849)
      #  elementName = 'Nacc'

    ###
    # Some Hydrogens are aromatic.
    # Based on probe:select.c:setAromaticProp()
    if elementName == 'H':
      if IsAromatic(resName, atomPadded):
        elementName = 'Har'

    ###
    # Some Hydrogens are polar
    # Based on probe:probe.c:updateHydrogenInfo()
    # If the Hydrogen is bonded to an N, O, or S, then mark it as polar and mark the
    # parent as not a donor.
    # @todo

    ###
    # Look up the element name, which fails if it is not in the table.
    try:
      ai = AtomInfo(self._Index[elementName])
    except Exception:
      return (None, "WARNING: atom "+atom.name+" from "+atom.parent().resname+
                 ' not recognized or unknown')

    # Return the value, warning if we've been asked to.
    warning = ""
    if emitWarning:
      warning = ("WARNING: atom "+atom.name+" from "+atom.parent().resname+
                 ' will be treated as '+elementName)
    return ( ai, warning )

  def _FindProperRadius(self, ai):
    """Given an AtomInfo in the atom table, find the appropriate radius.
    :param ai: AtomInfo to look up the radius in.
    :returns Proper radius, depending on construtor parameters.
    """
    if self._useNeutronDistances:
      return ai.vdwNeutronExplicit
    elif self._useImplicitHydrogenDistances:
      return ai.vdwElectronCloudImplicit
    else:
      return ai.vdwElectronCloudExplicit

  def FindProbeExtraAtomInfo(self, atom):
    """Given an iotbx.pdb.atom, look up its mmtbx_probe_ext.ExtraAtomInfo in the atom table.
    Note: Makes use of the mmtbx.probe.useNeutronDistances option to determine whether to
    return electron-cloud distance (default, when False) or neutron distances (when True).
    :param atom: iotbx.pdb.atom entry to look up.
    :returns a pair (mmtbx_probe_ext.ExtraAtomInfo structure filled with the info from the table,
    warning string) on success, raises ValueError on failure.  The warning string is empty
    if there is no warning, contains a printable warning if there was.
    """

    ai, warn = self.FindAtomInfo(atom)
    if ai is None:
      raise ValueError("FindProbeExtraAtomInfo(): Could not look up atom",atom.name)
    ret = probe.ExtraAtomInfo()
    ret.isAcceptor = bool(ai.flags & AtomFlags.ACCEPTOR_ATOM)
    ret.isDonor = bool(ai.flags & AtomFlags.DONOR_ATOM)
    ret.isDummyHydrogen = bool(ai.flags & AtomFlags.HB_ONLY_DUMMY_ATOM)
    ret.vdwRadius = self._FindProperRadius(ai)

    return ( ret, warn )

  def MaximumVDWRadius(self):
    """Return the maximum VdW radius of any atom type in our table.
    Cache the result after the first computation so that is faster when called
    more than once.
    """
    try:
      return self._maxVDW
    except Exception:
      self._maxVDW = self._FindProperRadius(AtomInfo(self._AtomTable[0]))
      for a in self._AtomTable[1:]:
        v = self._FindProperRadius(AtomInfo(a))
        if v > self._maxVDW:
          self._maxVDW = v
      return self._maxVDW

def Test(inFileName = None):

  #========================================================================
  # Make sure we can fill in mmtbx.probe.ExtraAtomInfoList info.
  # Generate an example data model with a small molecule in it unless we
  # were given a file name on the command line.
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    dm = DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
    mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()             #   get the model

  # Fill in an ExtraAtomInfoList with an entry for each atom in the hierarchy.
  # We first find the largest i_seq sequence number in the model and reserve that
  # many entries so we will always be able to fill in the entry for an atom.
  atoms = model.get_atoms()
  maxI = atoms[0].i_seq
  for a in atoms:
    if a.i_seq > maxI:
      maxI = a.i_seq
  extra = []
  for i in range(maxI+1):
    extra.append(probe.ExtraAtomInfo())

  # Traverse the hierarchy and look up the extra data to be filled in.
  # Get a list of all the atoms in the chain while we're at it
  at = AtomTypes()
  ph = model.get_hierarchy()
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          for a in ag.atoms():
            ei, warn = at.FindProbeExtraAtomInfo(a)
            extra[a.i_seq] = ei
            # User code should test for and print warnings
            #if len(warn) > 0:
            #  print(warn)

  #========================================================================
  # Find an Oxygen atom and ask for its radii with explicit Hydrogen, implicit Hydrogen,
  # and Nuclear radii.
  o = None
  ph = model.get_hierarchy()
  for a in ph.models()[0].atoms():
    if a.element.strip() == 'O':
      o = a
  assert o is not None, "AtomTypes.Test(): Could not find Oxygen (internal test failure)"
  explicitH = AtomTypes(useNeutronDistances = False,
                        useImplicitHydrogenDistances = False).FindProbeExtraAtomInfo(o)[0].vdwRadius
  implicitH = AtomTypes(useNeutronDistances = False,
                        useImplicitHydrogenDistances = True).FindProbeExtraAtomInfo(o)[0].vdwRadius
  neutronH = AtomTypes(useNeutronDistances = True,
                        useImplicitHydrogenDistances = False).FindProbeExtraAtomInfo(o)[0].vdwRadius
  assert explicitH != implicitH, "AtomTypes.Test(): Implicit and explicit Oxygen radii did not differ as expected"

  #========================================================================
  # Check MaximumVDWRadius, calling it twice to make sure both the cached and non-cached
  # results work.
  for i in range(2):
    assert at.MaximumVDWRadius() == 2.5, "AtomTypes.Test(): Unexpected MaximumVDWRadius(): got "+str(MaximumVDWRadius())+", expected 2.5"

  #========================================================================
  # Check IsAromatic() to ensure it gives results when expected and not when not.
  aromaticChecks = [
      ['PHE', 'CE2', True],
      ['  U', 'HN3', True],
      ['ASN', 'O', False]
    ]
  for a in aromaticChecks:
    assert IsAromatic(a[0],a[1]) == a[2], "AtomTypes.Test(): {} {} not marked as aromatic {}".format(a[0],a[1],a[2])

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
  print('OK')
