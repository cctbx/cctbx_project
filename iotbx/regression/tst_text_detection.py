from iotbx.file_io.text_detection import (
  sniff_text_datatype, _looks_like_pdb, _looks_like_ncs_spec,
  _looks_like_phil, _looks_like_sequence)

PDB = (
  'CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1\n'
  'ATOM      1  N   GLY A   1      -9.0   4.6   6.1  1.00 16.8           N\n'
  'END\n')

NCS = (
  '\nSummary of NCS information\nThu May 21 14:48:36 2020\n\n\n'
  'new_ncs_group\nnew_operator\n\n'
  'rota_matrix    1.0000    0.0000    0.0000\n'
  'rota_matrix    0.0000    1.0000    0.0000\n'
  'rota_matrix    0.0000    0.0000    1.0000\n'
  'tran_orth     0.0000    0.0000    0.0000\n\n'
  'center_orth   14.4035    7.4690    0.2556\n'
  'CHAIN A\nRMSD 0\nMATCHING 3\n  RESSEQ 1:3\n')

PHIL_NCS = (
  'refinement.pdb_interpretation.ncs_group {\n'
  '  reference = chain A and resseq 20:28\n'
  '  selection = chain B and resseq 20:28\n}\n')

PHIL_SCOPE_ONLY = 'refinement.pdb_interpretation.ncs_group {\n}\n'
PHIL_BASIC = 'foo {\n  bar = 1\n    .type = int\n}\n'
PHIL_COMMENTED = '# a leading comment\n# another\nmaster.thing = 5\n'
SEQ_FASTA = '>1SAR\nDVSGTVCLSALPPEATDTLNLIASDGPFPYSQDG\n'
SEQ_BARE = 'DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDG\n'
STRAY_RECORDS = 'foo = 1\nHEADER text at column zero\nSHEET more text\n'
PROSE = 'This is explanatory prose.\nNothing structured here.\n'
DAT_NUMERIC = '1.0 2.0 3.0\n4.0 5.0 6.0\n'
HEADER_ONLY = (                       # crystal symmetry but no atoms -> not model
  'CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1\n'
  'SCALE1      0.045585  0.000000  0.014006        0.00000\n')
GAPPED_SEQ = '>aln\nDVSG--TVCLS--ALPP\n'   # alignment gaps -> any_file reads 'aln'
LONE_HEADER = '>just a header, no residues\n'

def test_pdb_predicate():
  assert _looks_like_pdb(PDB) is True
  assert _looks_like_pdb(PHIL_NCS) is False
  assert _looks_like_pdb(SEQ_BARE) is False
  # word-like stray records do NOT count as model
  assert _looks_like_pdb(STRAY_RECORDS) is False
  # an atom-bearing record is required: a header/symmetry-only file is not a
  # model (any_file._try_as_pdb needs models_size() > 0)
  assert _looks_like_pdb(HEADER_ONLY) is False
  print('test_pdb_predicate OK')

def test_ncs_spec_predicate():
  assert _looks_like_ncs_spec(NCS) is True
  # a PHIL ncs_group scope is NOT an ncs_spec (we key on new_ncs_group)
  assert _looks_like_ncs_spec(PHIL_NCS) is False
  assert _looks_like_ncs_spec(PDB) is False
  print('test_ncs_spec_predicate OK')

def test_phil_predicate():
  assert _looks_like_phil(PHIL_NCS) is True
  assert _looks_like_phil(PHIL_BASIC) is True
  # dotted scope-open with no '=' lines must still match
  assert _looks_like_phil(PHIL_SCOPE_ONLY) is True
  assert _looks_like_phil(PHIL_COMMENTED) is True
  # an ncs_spec body and a residue block are NOT phil
  assert _looks_like_phil(NCS) is False
  assert _looks_like_phil(SEQ_BARE) is False
  print('test_phil_predicate OK')

def test_sequence_predicate():
  assert _looks_like_sequence(SEQ_FASTA) is True
  assert _looks_like_sequence(SEQ_BARE) is True
  assert _looks_like_sequence(DAT_NUMERIC) is False
  assert _looks_like_sequence(PHIL_BASIC) is False
  # gaps and lone headers do NOT qualify (any_file reads these as 'aln', not seq)
  assert _looks_like_sequence(GAPPED_SEQ) is False
  assert _looks_like_sequence(LONE_HEADER) is False
  print('test_sequence_predicate OK')

def test_precedence_and_reclassify():
  assert sniff_text_datatype(PDB) == 'model'
  assert sniff_text_datatype(NCS) == 'ncs_spec'
  # the .ncs-holds-PHIL case: content wins -> phil (not ncs_spec)
  assert sniff_text_datatype(PHIL_NCS) == 'phil'
  assert sniff_text_datatype(PHIL_SCOPE_ONLY) == 'phil'
  assert sniff_text_datatype(SEQ_FASTA) == 'sequence'
  # a stray HEADER/SHEET line does not win 'model'; the '=' line wins 'phil'
  assert sniff_text_datatype(STRAY_RECORDS) == 'phil'
  # header/symmetry-only PDB and a gapped sequence are not reclassified
  assert sniff_text_datatype(HEADER_ONLY) is None
  assert sniff_text_datatype(GAPPED_SEQ) is None
  print('test_precedence_and_reclassify OK')

def test_no_match():
  assert sniff_text_datatype(PROSE) is None
  assert sniff_text_datatype(DAT_NUMERIC) is None
  assert sniff_text_datatype('') is None
  print('test_no_match OK')

if __name__ == '__main__':
  test_pdb_predicate()
  test_ncs_spec_predicate()
  test_phil_predicate()
  test_sequence_predicate()
  test_precedence_and_reclassify()
  test_no_match()
  print('OK')
