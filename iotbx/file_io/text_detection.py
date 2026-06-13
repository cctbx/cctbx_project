'''
iotbx.file_io.text_detection: positive content verification for the text
DataManager datatypes whose extensions collide or bind weakly -- phil, ncs_spec,
model (PDB format), and sequence.

sniff_text_datatype(text) classifies decoded text by cheap, anchored structural
signatures: the text counterpart to the binary magic checks in detection.py. It
never parses (read_file remains the single parse authority) and never raises on
str input. detection.get_file_type calls it so text detection is
content-authoritative among the sniffable types -- a confident content match
wins over the extension.
'''
import re

# -- model (PDB format) -------------------------------------------------------
# Require an atom-bearing coordinate record (ATOM/HETATM) anchored at column 0.
# any_file._try_as_pdb only accepts a file with models_size() > 0, so a
# header-only file (CRYST1/SCALE/ORIGX with no atoms) is NOT a model to the
# reader; detection must agree, or it reclassifies a non-model and the following
# read fails. Header/symmetry records alone are deliberately insufficient, and
# word-like records (HEADER/HELIX/SHEET/...) never qualify. (mmCIF models arrive
# as .cif -> CIF token-scan, not here.)
_PDB_RECORD = re.compile(r'^(ATOM  |HETATM)')

def _looks_like_pdb(text):
  '''
  Whether text carries an atom-bearing column-0 PDB record (ATOM/HETATM).

  Parameters
  ----------
  text : str
      Decoded file prefix.

  Returns
  -------
  bool
      True if any line begins with an ATOM or HETATM record.
  '''
  for line in text.splitlines():
    if _PDB_RECORD.match(line):
      return True
  return False

# -- ncs_spec -----------------------------------------------------------------
# Structural tokens emitted by mmtbx.ncs.ncs.as_ncs_spec_string. Distinct from a
# PHIL `ncs_group {` scope: we key on new_ncs_group, not bare ncs_group, so a
# PHIL ncs_group scope does not match here.
_NCS_RE = re.compile(
  r'\b(new_ncs_group|new_operator|rota_matrix|tran_orth|center_orth)\b')

def _looks_like_ncs_spec(text):
  '''
  Whether text carries ncs_spec structural tokens.

  Parameters
  ----------
  text : str
      Decoded file prefix.

  Returns
  -------
  bool
      True if a new_ncs_group / new_operator / rota_matrix / tran_orth /
      center_orth token is present.
  '''
  return _NCS_RE.search(text) is not None

# -- phil ---------------------------------------------------------------------
# A phil definition (dotted name = value) or a scope open (dotted name {). The
# dotted prefix is required so refinement.pdb_interpretation.ncs_group { matches.
_PHIL_DEF = re.compile(r'^\s*[A-Za-z_]\w*(\.[A-Za-z_]\w*)*\s*=')
_PHIL_SCOPE = re.compile(r'^\s*[A-Za-z_]\w*(\.[A-Za-z_]\w*)*\s*\{')

def _looks_like_phil(text):
  '''
  Whether text shows phil structure, ignoring full-line # comments.

  Parameters
  ----------
  text : str
      Decoded file prefix.

  Returns
  -------
  bool
      True if a non-comment line is a phil definition (dotted name = value) or a
      scope open (dotted name {).
  '''
  for line in text.splitlines():
    stripped = line.strip()
    if (not stripped) or stripped.startswith('#'):
      continue
    if _PHIL_DEF.match(line) or _PHIL_SCOPE.match(line):
      return True
  return False

# -- sequence -----------------------------------------------------------------
# Residue-letter content, optionally with FASTA/PIR '>' header lines. Gaps ('-')
# are excluded and a header alone (no residues) does not qualify: any_file's
# _try_as_seq reads gapped/alignment input as 'aln' (not 'seq') and needs actual
# residues, so detection matches it. Loosest signal; runs last so it only ever
# claims text the other predicates declined.
_RESIDUE = set('ACDEFGHIKLMNPQRSTVWYBXZU*')

def _looks_like_sequence(text):
  '''
  Whether text is residue-letter content (no gaps), ignoring FASTA headers.

  Parameters
  ----------
  text : str
      Decoded file prefix.

  Returns
  -------
  bool
      True if every non-header line is IUPAC residue letters (no gaps), with at
      least one such line.
  '''
  nonblank = [line.strip() for line in text.splitlines() if line.strip()]
  if not nonblank:
    return False
  residue_lines = [line for line in nonblank if not line.startswith('>')]
  if not residue_lines:           # only header(s), no residue content
    return False
  for line in residue_lines:
    for ch in line:
      if (ch.upper() not in _RESIDUE) and (not ch.isspace()):
        return False
  return True

# -- entry point --------------------------------------------------------------
_PRECEDENCE = (
  ('model',    _looks_like_pdb),
  ('ncs_spec', _looks_like_ncs_spec),
  ('phil',     _looks_like_phil),
  ('sequence', _looks_like_sequence),
)

def sniff_text_datatype(text):
  '''
  Classify decoded text content by structural signature.

  Parameters
  ----------
  text : str
      Decoded file prefix.

  Returns
  -------
  str or None
      'model', 'ncs_spec', 'phil', or 'sequence' for the highest-precedence
      signature that matches; None if none match. A pure content classifier --
      valid_types filtering is the caller's (get_file_type's) job.
  '''
  for datatype, predicate in _PRECEDENCE:
    if predicate(text):
      return datatype
  return None
