from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.utils import Sorry
from mmtbx.geometry_restraints.torsion_restraints.utils import \
    check_for_internal_chain_ter_records

# check_for_internal_chain_ter_records() guards the torsion / reference-model
# machinery against TER cards that split a single contiguous chain in two.
# It must NOT flag a TER that legitimately terminates a chain whose id is
# reused later in the file -- a very common layout where ligands (or a second
# discontinuous segment) sharing a protein chain id are grouped at the end of
# the file (e.g. PDB 1rhq, whose peptide-like 0ZZ ligand reuses chain ids
# 'A'/'D' and is written after all four protein chains).

def _prot_chain(chain_id, resseq1, resseq2, x0):
  # Two amino-acid residues (GLY, ALA) forming a minimal protein chain.
  lines = [
    "ATOM      1  N   GLY %s%4d    %8.3f  11.000  11.000  1.00 20.00           N" % (chain_id, resseq1, x0),
    "ATOM      2  CA  GLY %s%4d    %8.3f  12.000  11.000  1.00 20.00           C" % (chain_id, resseq1, x0+1),
    "ATOM      3  C   GLY %s%4d    %8.3f  13.000  11.000  1.00 20.00           C" % (chain_id, resseq1, x0+2),
    "ATOM      4  O   GLY %s%4d    %8.3f  14.000  11.000  1.00 20.00           O" % (chain_id, resseq1, x0+3),
    "ATOM      5  N   ALA %s%4d    %8.3f  15.000  11.000  1.00 20.00           N" % (chain_id, resseq2, x0),
    "ATOM      6  CA  ALA %s%4d    %8.3f  16.000  11.000  1.00 20.00           C" % (chain_id, resseq2, x0+1),
    "ATOM      7  C   ALA %s%4d    %8.3f  17.000  11.000  1.00 20.00           C" % (chain_id, resseq2, x0+2),
    "ATOM      8  O   ALA %s%4d    %8.3f  18.000  11.000  1.00 20.00           O" % (chain_id, resseq2, x0+3),
    ]
  return "\n".join(lines) + "\n"

CRYST = "CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1\n"

def _raises_errant_ter(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str.splitlines())
  h = pdb_inp.construct_hierarchy()
  h.reset_i_seq_if_necessary()
  try:
    check_for_internal_chain_ter_records(
        pdb_hierarchy=h, ter_indices=pdb_inp.ter_indices())
  except Sorry:
    return True
  return False

def exercise_reused_chain_id_at_end():
  # Protein chain 'A', legitimate TER, protein chain 'B', legitimate TER, then
  # a second protein segment that reuses chain id 'A' at the end of the file.
  # The two TER cards both terminate a chain properly; none is internal.
  pdb_str = (CRYST
      + _prot_chain("A", 1, 2, 11.0) + "TER\n"
      + _prot_chain("B", 1, 2, 31.0) + "TER\n"
      + _prot_chain("A", 10, 11, 51.0) + "END\n")
  assert not _raises_errant_ter(pdb_str), \
      "legitimate TER before a reused chain id was flagged as errant"

def exercise_genuine_internal_ter():
  # A single contiguous protein chain 'A' broken in the middle by a spurious
  # TER. This must still be detected.
  pdb_str = (CRYST
      + "ATOM      1  N   GLY A   1      11.000  11.000  11.000  1.00 20.00           N\n"
      + "ATOM      2  CA  GLY A   1      12.000  11.000  11.000  1.00 20.00           C\n"
      + "ATOM      3  C   GLY A   1      13.000  11.000  11.000  1.00 20.00           C\n"
      + "ATOM      4  O   GLY A   1      14.000  11.000  11.000  1.00 20.00           O\n"
      + "TER\n"
      + "ATOM      5  N   ALA A   2      15.000  11.000  11.000  1.00 20.00           N\n"
      + "ATOM      6  CA  ALA A   2      16.000  11.000  11.000  1.00 20.00           C\n"
      + "ATOM      7  C   ALA A   2      17.000  11.000  11.000  1.00 20.00           C\n"
      + "ATOM      8  O   ALA A   2      18.000  11.000  11.000  1.00 20.00           O\n"
      + "END\n")
  assert _raises_errant_ter(pdb_str), \
      "a spurious TER inside a single contiguous chain was not detected"

def exercise_contiguous_non_protein_ligand():
  # Protein chain 'A', TER, then a non-protein ligand that reuses chain id 'A'
  # immediately afterwards. The TER terminates the protein legitimately.
  pdb_str = (CRYST
      + _prot_chain("A", 1, 2, 11.0) + "TER\n"
      + "HETATM   19  C1  EDO A 101      30.000  30.000  30.000  1.00 20.00           C\n"
      + "HETATM   20  O1  EDO A 101      31.000  30.000  30.000  1.00 20.00           O\n"
      + "END\n")
  assert not _raises_errant_ter(pdb_str), \
      "TER before a contiguous non-protein ligand was flagged as errant"

def exercise_no_ter():
  pdb_str = CRYST + _prot_chain("A", 1, 2, 11.0) + "END\n"
  assert not _raises_errant_ter(pdb_str)

def run():
  exercise_reused_chain_id_at_end()
  exercise_genuine_internal_ter()
  exercise_contiguous_non_protein_ligand()
  exercise_no_ter()
  print("OK")

if (__name__ == "__main__"):
  run()
