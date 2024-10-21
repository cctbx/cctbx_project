from __future__ import absolute_import, division, print_function
from iotbx.pdb import amino_acid_codes, nucleic_acid_codes
from scitbx.matrix import col

_amino_acid_resnames = sorted(amino_acid_codes.one_letter_given_three_letter.keys())
def _IsStandardResidue(resname):
  return resname.strip().upper() in _amino_acid_resnames

# Find the RNA and DNA residue sets. Remove the RNA names from the DNA set to get only definitely DNA names.
_rna_resnames = set(nucleic_acid_codes.rna_one_letter_code_dict.keys())
_dna_resnames = set(nucleic_acid_codes.dna_one_letter_code_dict.keys()) - _rna_resnames
_nucleic_acid_resnames = _dna_resnames.union(_rna_resnames)
def _IsNucleicAcidResidue(resname):
  return resname.strip().upper() in _nucleic_acid_resnames

def _FindNamedAtomInResidue(residue_group, atom_names):
  '''Return the first atom in the residue with a name in the atom_names list.
  :param residue_group: iotbx.pdb.hierarchy.residue_group, the residue to search.
  :param atom_names: list of strings, the names of the atoms to search for.
  :return: iotbx.pdb.hierarchy.atom, the first atom found with a name in the atom_names list,
  or None if none are found.
  '''
  for name in atom_names:
    name = name.strip().upper()
    for atom in residue_group.atoms():
      if atom.name.strip().upper() == name.upper():
        return atom
  return None

def _FindContiguousResiduesByAtomDistances(chain, type_function, desired_atoms, distance_threshold):
  '''Return a list of contiguous nucleic acid residues in the chain based on the distance between P atoms,
  or O5* or O5' if the P atom is not found.
  :param chain: PDB chain to be searched for contiguous residues.
  :param type_function: function, a function that takes a residue name and returns True if the residue is of the desired type.
  :param desired_atoms: list of strings, the names of the atoms to use for distance calculations.  The first
  one found in the residue will be used.
  :param distance_threshold: float, the maximum distance between desired atoms to consider them contiguous.
  '''
  contiguous_residues = []

  current_contig_residues = []
  prev_atom = None

  for residue_group in chain.residue_groups():
    # Check if the residue is of the desired type, based on its name
    if not type_function(residue_group.unique_resnames()[0]):
      continue

    # Attempt to find the desired atom in the current residue.  Search all atoms for the
    # first one in the desired_atoms list, then later ones if it is not found.
    # If no atom is found, skip to the next residue
    atom = _FindNamedAtomInResidue(residue_group, desired_atoms)
    if atom is None:
      continue

    # If this is the first residue being examined, initialize the list
    if prev_atom is None:
      current_contig_residues.append(residue_group)
    else:
      # Calculate the distance between the current and previous atoms
      current_pos = col(atom.xyz)
      prev_pos = col(prev_atom.xyz)
      distance = (current_pos - prev_pos).length()

      if distance < distance_threshold:
        # If the distance is within the threshold, add to the current list
        current_contig_residues.append(residue_group)
      else:
        # If not, start a new list for the current residue if we have at least one residue
        if len(current_contig_residues) > 0:
          contiguous_residues.append(current_contig_residues)
        current_contig_residues = [residue_group]

    # Update the previous atom to the current one
    prev_atom = atom

  # After iterating through the chain, add any remaining contiguous residues to the main list
  # if there are at least two.
  if len(current_contig_residues) > 1:
    contiguous_residues.append(current_contig_residues)

  return contiguous_residues

# ------------------------------------------------------------------------------

def chain_has_DNA(chain):
  '''Return True if the chain contains any DNA residues.
  :param chain: PDB chain to be searched for DNA residues.
  '''
  for residue_group in chain.residue_groups():
    if residue_group.unique_resnames()[0].strip().upper() in _dna_resnames:
      return True
  return False

def chain_has_RNA(chain):
  '''Return True if the chain contains any RNA residues.
  :param chain: PDB chain to be searched for RNA residues.
  '''
  for residue_group in chain.residue_groups():
    if residue_group.unique_resnames()[0].strip().upper() in _rna_resnames:
      return True
  return False

# ------------------------------------------------------------------------------

def find_contiguous_protein_residues(chain, distance_threshold=5.0):
  '''Return a list of contiguous protein residues in the chain based on the distance between CA atoms.
  :param chain: PDB chain to be searched for contiguous residues.
  :param distance_threshold: float, the maximum distance between CA atoms to consider them contiguous.
  The empirical 5.0 default value comes from the Richardson lab's Prekin code; the ideal length is 3.80.
  '''
  return _FindContiguousResiduesByAtomDistances(chain, _IsStandardResidue, ["CA"], distance_threshold)

def find_contiguous_nucleic_acid_residues(chain, distance_threshold=10.0):
  '''Return a list of contiguous nucleic acid residues in the chain based on the distance between P atoms,
  or O5* or O5' if the P atom is not found.
  :param chain: PDB chain to be searched for contiguous residues.
  :param distance_threshold: float, the maximum distance between CA atoms to consider them contiguous.
  The empirical 10.0 default value comes from the Richardson lab's Prekin code; the ideal length is ~7.
  '''
  return _FindContiguousResiduesByAtomDistances(chain, _IsNucleicAcidResidue, ["P", "O5*", "O5'"], distance_threshold)

# ------------------------------------------------------------------------------

class GuidePoint:
  '''Class to hold a guide point for the ribbons representation.'''
  def __init__(self, pos = None, cvec = None, dvec = None, offsetFactor = 0.0, widthFactor = 0.0, prevRes = None, nextRes = None):
    '''Constructor for the class.
    :param pos: vec3_double, the final position of the guide point.
    :param cvec: vec3_double, the unit vector normal to the local plane of the ribbon.
    :param dvec: vec3_double, the unit vector in the plane of the ribbon, perpendicular to its overall direction.
    :param offsetFactor: float, the factor needed to pull spline through guidepoints.
    It is already applied and does not depend on neighbors (unlike widthFactor).
    :param widthFactor: float, Suggested ribbon width modifier, from 0 (skinniest) to 1 (fattest).
    This is based on low/high curvatire of 3+ residues in a row.
    :param prevRes: iotbx.pdb.hierarchy.residue_group, the residue in the chain just before this guidepoint
    (never None). Usually different residues for proteins and for nucleic acids, except at chain ends.
    For nucleic acids, prevRes is the residue centered on the guide point and nextRes is the one after.
    (Protein guides fall between residues; nucleic guides fall in middle of a residue.)
    :param nextRes: iotbx.pdb.hierarchy.residue_group, the residue in the chain just after this guidepoint.
    '''
    self.pos = pos
    self.cvec = cvec
    self.dvec = dvec
    self.offsetFactor = offsetFactor
    self.widthFactor = widthFactor
    self.prevRes = prevRes
    self.nextRes = nextRes

  def __str__(self):
    return 'GuidePoint(pos={}, cvec={}, dvec={}, offsetFactor={}, widthFactor={}, prev={}, next={})'.format(
      self.pos, self.cvec, self.dvec, self.offsetFactor, self.widthFactor, self.prevRes, self.nextRes)

# ------------------------------------------------------------------------------

def make_protein_guidepoints(contiguous_residues):
  '''Return a list of GuidePoint objects for the protein residues.
  :param contiguous_residues: list of iotbx.pdb.hierarchy.residue_group entries, on of
  the contiguous residue lists returned by find_contiguous_protein_residues().
  '''

  # Initialize an empty list of guiepoints that lie between residues (one less than the number of residues)
  # plus two at the beginning and two at the end.  These will be filled in below.  These mist be independent
  # GuidePoint objects, not references to the same object, because they will be modified in place.
  numResidues = len(contiguous_residues)
  guidepoints = []
  for i in range(numResidues - 1 + 2 + 2):
    guidepoints.append(GuidePoint())

  # Initialize some values for the calculations
  maxOffset = 1.5   # Maximum displacement of guidepoint based on curvature
  anchorAtoms = ["CA"]  # Atoms to use for anchor points

  # Make normal guidepoints between each pair of residues.
  for i in range(numResidues - 1):

    # Find the guidepoint to work on.
    g = guidepoints[i+2]

    # Set the residues before and after the guidepoint
    g.prevRes = contiguous_residues[i]
    g.nextRes = contiguous_residues[i+1]

    # Find the anchor atoms for the current and previous residues.
    # We know they are there because we found the residues in the contiguous_residues list.
    ca1 = _FindNamedAtomInResidue(g.prevRes, anchorAtoms)
    ca2 = _FindNamedAtomInResidue(g.nextRes, anchorAtoms)

    # Calculate the position of the guide point, which is halfway between them.
    g.pos = (col(ca1.xyz) + col(ca2.xyz)) / 2

    # Based on Ca(i-1) to Ca(i+2) distance, we may adjust ribbon width
    # and/or guidepoint position. Ribbon widens in areas of high OR low
    # curvature (alpha/beta, respectively). This is only a preliminary
    # estimate -- it's applied only for 3+ residues in a row. (checked below)
    #
    # For high curvature ONLY (alpha), we offset the guidepoint
    # outwards to make the spline track the chain (esp. helix);
    # this is done for each and does not require 3+ in a row.
    # Offset vector goes from the midpoint of Ca(i-1) to Ca(i+2)
    # thru the current guide point
    # (which is the midpoint of Ca(i) and Ca(i+1)).
    #
    # CA-CA DIST  WIDTH FACTOR    OFFSET FACTOR   NOTE
    # ==========  ============    =============   ====
    # 5.0         1               1               ~limit for curled-up protein
    # 5.5         1               1               } linear interpolation
    # 7.0         0               0               } from 1.0 to 0.0
    # 9.0         0               0               } linear interpolation
    # 10.5        1               0               } from 1.0 to 0.0
    # 11.0        1               0               ~limit for extended protein

    if 1 <= i and i <= numResidues - 3:
      ca0 = _FindNamedAtomInResidue(contiguous_residues[i-1], anchorAtoms)
      ca3 = _FindNamedAtomInResidue(contiguous_residues[i+2], anchorAtoms)
      cacaDist = (col(ca0.xyz) - col(ca3.xyz)).length()
      if cacaDist < 7:
        g.widthFactor = min(1.5, 7-cacaDist) / 1.5
        g.offsetFactor = g.widthFactor
        midpt = (col(ca0.xyz) + col(ca3.xyz)) / 2
        offsetVec = (g.pos - midpt).normalize()
        g.pos += maxOffset * g.offsetFactor * offsetVec
      elif cacaDist > 9:
        g.widthFactor = min(1.5, cacaDist-9) / 1.5
        g.offsetFactor = 0
      else:
        g.widthFactor = 0
        g.offsetFactor = 0

    # Do this last so that for CA-only structures, everything above has been calculated
    # before we might have to bail on this guidepoint.
    ox1 = _FindNamedAtomInResidue(g.prevRes, ["O"])
    if ox1 is not None:
      avec = col(ca2.xyz) - col(ca1.xyz)
      bvec = col(ox1.xyz) - col(ca1.xyz)
      g.cvec = avec.cross(bvec).normalize()
      g.dvec = g.cvec.cross(avec).normalize()
    else:
      g.cvec = col({0, 0, 0})
      g.dvec = col({0, 0, 0})

  # Check on widthFactors -- only apply for 3+ in a row > 0
  i = 2
  while i < len(guidepoints) - 2:
    # Scan to find the first widened guidepoint.
    if guidepoints[i].widthFactor == 0:
      i += 1
      continue
    # Scan to find the last widened guidepoint.
    firstWide = i
    nextThin = i+1
    while nextThin < len(guidepoints) - 2 and guidepoints[nextThin].widthFactor != 0:
      nextThin += 1
    # If the span is less than 3, set them all back to zero.
    if nextThin - firstWide < 3:
      for j in range(firstWide, nextThin):
        guidepoints[j].widthFactor = 0
    i += 1

  # Make two dummy guidepoints at the beginning and end of the chain.
  firstGuides = GuidePoint(col(_FindNamedAtomInResidue(contiguous_residues[0], anchorAtoms).xyz),
                           guidepoints[2].cvec, guidepoints[2].dvec,
                           guidepoints[2].offsetFactor, guidepoints[2].widthFactor,
                           contiguous_residues[0], contiguous_residues[0])
  guidepoints[0] = firstGuides
  guidepoints[1] = firstGuides
  lastGuides = GuidePoint(col(_FindNamedAtomInResidue(contiguous_residues[-1], anchorAtoms).xyz),
                          guidepoints[-3].cvec, guidepoints[-3].dvec,
                          guidepoints[-3].offsetFactor, guidepoints[-3].widthFactor,
                          contiguous_residues[-1], contiguous_residues[-1])
  guidepoints[-1] = lastGuides
  guidepoints[-2] = lastGuides

  return guidepoints

def make_nucleic_acid_guidepoints(contiguous_residues):
  '''Return a list of GuidePoint objects for the nucleic acid residues.
  :param contiguous_residues: list of iotbx.pdb.hierarchy.residue_group entries, one of
 the contiguous residue lists returned by find_contiguous_amino_acid_residues().
  '''

  # Initialize an empty list of guiepoints that lie between residues (one less than the number of residues)
  # plus two at the beginning and two at the end.  These will be filled in below.  These mist be independent
  # GuidePoint objects, not references to the same object, because they will be modified in place.
  numResidues = len(contiguous_residues)
  guidepoints = []
  for i in range(numResidues - 1 + 2 + 2):
    guidepoints.append(GuidePoint())

  # Initialize some values for the calculations
  maxOffset = 1.5   # Maximum displacement of guidepoint based on curvature
  anchorAtoms = ["P", "O5*", "O5'"]  # Atoms to use for anchor points

  # Make normal guidepoints between each pair of residues.
  for i in range(numResidues - 1):

    # Find the guidepoint to work on.
    g = guidepoints[i+2]

    # Set the residues before and after the guidepoint
    g.prevRes = contiguous_residues[i]
    g.nextRes = contiguous_residues[i+1]

    # Find the anchor atoms for the current and previous residues.
    # We know they are there because we found the residues in the contiguous_residues list.
    phos1 = _FindNamedAtomInResidue(g.prevRes, anchorAtoms)
    phos2 = _FindNamedAtomInResidue(g.nextRes, anchorAtoms)

    # Calculate the position of the guide point, which is halfway between them.
    g.pos = (col(phos1.xyz) + col(phos2.xyz)) / 2

    # Based on P(i-1) to P(i+1) or P(i) to P(i+2) distance, we may offset guide points.
    # For areas of high curvature, the guide point is offset towards the C4'Based
    # on P(i-1) to P(i+2) distance, we may adjust ribbon width.

    if (1 <= i and i <= numResidues - 3):
      p0 = _FindNamedAtomInResidue(contiguous_residues[i-1], anchorAtoms)
      p3 = _FindNamedAtomInResidue(contiguous_residues[i+2], anchorAtoms)
      ppDist1 = (g.pos - col(p0.xyz)).length()
      ppDist2 = (g.pos - col(p3.xyz)).length()

      # Default values for regions of low curvature
      isTightlyCurved = False
      ppDistance = ppDist1
      maxOffset = 0.0
      g.offsetFactor = 0.0
      g.widthFactor = 0.0   # There is never width hinting for nucleic acids
      if ppDist1 <= 9.0 or ppDist2 <= 9.0:
        isTightlyCurved = True
        ppDistance = min(ppDist1, ppDist2)
        carbon4 = _FindNamedAtomInResidue(g.prevRes, ["C4*", "C4'"])
        if carbon4 is not None:
          maxOffset = (g.pos - col(carbon4.xyz)).length() + 1.0  # Allows guidepoint to go past the C4'
          g.offsetFactor = (9.0 - ppDistance) / (9.0 - 7.0)
          g.offsetFactor = min(1.0, g.offsetFactor)   # Reaches full offset at 7A curvature
          offsetVec = (col(carbon4.xyz) - g.pos).normalize()
          g.pos += maxOffset * g.offsetFactor * offsetVec

    # We do this last so that if there are missing atoms everying above has been calculated
    # before we might have to bail on this guidepoint.
    carbon3 = _FindNamedAtomInResidue(g.prevRes, ["C3*", "C3'"])
    carbon1 = _FindNamedAtomInResidue(g.nextRes, ["C1*", "C1'"])
    if carbon3 is not None and carbon1 is not None:
      avec = col(phos2.xyz) - col(phos1.xyz)
      bvec = col(carbon1.xyz) - col(carbon3.xyz)
      g.cvec = avec.cross(bvec).normalize()
      g.dvec = g.cvec.cross(avec).normalize()
    else:
      g.cvec = col({0, 0, 0})
      g.dvec = col({0, 0, 0})

  # Make dummy guidepoints at the beginning and end of the chain.
  firstGuides = GuidePoint(col(_FindNamedAtomInResidue(contiguous_residues[0], anchorAtoms).xyz),
                            guidepoints[2].cvec, guidepoints[2].dvec,
                            guidepoints[2].offsetFactor, guidepoints[2].widthFactor,
                            contiguous_residues[0], contiguous_residues[0])
  guidepoints[0] = firstGuides
  guidepoints[1] = firstGuides
  # Prekin: 3' guide point is 2/3 of the way to the O3' from the last guide point,
  # but we just use the C3 atom.
  carbon3 = _FindNamedAtomInResidue(contiguous_residues[-1], ["C3*", "C3'", "P", "O5'", "O5*"])
  if carbon3 is not None:
    lastGuides = GuidePoint(col(carbon3.xyz),
                            guidepoints[-3].cvec, guidepoints[-3].dvec,
                            guidepoints[-3].offsetFactor, guidepoints[-3].widthFactor,
                            contiguous_residues[-1], contiguous_residues[-1])
    guidepoints[-1] = lastGuides
    guidepoints[-2] = lastGuides

  return guidepoints

# ------------------------------------------------------------------------------

def swap_edge_and_face(guidepoints):
  '''Converts a RNA-style ribbon (the default) to a DNA-style one,
    and vice versa, by swapping the "c" and "d" vectors.
    This results in a 90 degree rotation of the ribbon around the ribbon axis.
  :param guidepoints: list of GuidePoint objects to swap the edge and face vectors for.
  '''
  for g in guidepoints:
    g.cvec, g.dvec = g.dvec, g.cvec

def untwist_ribbon(guidepoints):
  '''Removes excess twist from a ribbon by reversing the sign of the GuidePoint
    "d" vectors as necessary. A vector is flipped if it has a negative dot product
    with the previous one (i.e. the angle between them is greater than 90 degrees).
  :param guidepoints: list of GuidePoint objects to untwist.
  '''
  for i in range(len(guidepoints)-1):
    if guidepoints[i].dvec.dot(guidepoints[i+1].dvec) < 0:
      guidepoints[i+1].dvec *= -1.0

# ------------------------------------------------------------------------------

def non_CA_atoms_present(structure):
  '''Return True if any atom other than CA is present in the structure.
  :param structure: iotbx.pdb.hierarchy.root object holding the structure to check.
  '''
  for model in structure.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom in residue_group.atoms():
          if atom.name.strip().upper() != "CA":
            return True
  return False
