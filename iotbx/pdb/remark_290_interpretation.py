from cctbx import sgtbx

example = """\
REMARK 290
REMARK 290
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY
REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 21 21 21
REMARK 290
REMARK 290      SYMOP   SYMMETRY
REMARK 290     NNNMMM   OPERATOR
REMARK 290       1555   X,Y,Z
REMARK 290       2555   1/2-X,-Y,1/2+Z
REMARK 290       3555   -X,1/2+Y,1/2-Z
REMARK 290       4555   1/2+X,1/2-Y,-Z
REMARK 290
REMARK 290     WHERE NNN -> OPERATOR NUMBER
REMARK 290           MMM -> TRANSLATION VECTOR
REMARK 290
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS
REMARK 290 THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM
REMARK 290 RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY
REMARK 290 RELATED MOLECULES.
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       36.30027
REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       59.50256
REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000
REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       46.45545
REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       59.50256
REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       36.30027
REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       46.45545
REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000
REMARK 290
REMARK 290 REMARK: NULL
"""

def _nnnmmm_operator(record, expected_nnn):
  if (not record.startswith("REMARK 290     ")): return None
  flds = record.split()
  if (len(flds) != 4): return None
  nnnmmm = flds[2]
  if (len(nnnmmm) < 4): return None
  nnn,mmm = nnnmmm[:-3], nnnmmm[-3:]
  if (mmm != "555"): return None
  try: nnn = int(nnn)
  except ValueError: return None
  if (nnn != expected_nnn): return None
  try: return sgtbx.rt_mx(flds[3])
  except RuntimeError: return None

def extract_symmetry_operators(remark_290_records):
  symmetry_operators = []
  status = 0
  for record in remark_290_records:
    if (record.startswith("REMARK 290     ")):
      symmetry_operator = _nnnmmm_operator(
        record=record,
        expected_nnn=len(symmetry_operators)+1)
      if (symmetry_operator is not None):
        if (status == 2): return None
        symmetry_operators.append(symmetry_operator)
        status = 1
      elif (status == 1):
        status == 2
  if (len(symmetry_operators) == 0): return None
  return symmetry_operators

def get_link_symmetry_operator(symmetry_operators, link_sym):
  if (symmetry_operators is None): return None
  link_sym = link_sym.strip()
  if (len(link_sym) < 4): return None
  if (link_sym[-4] == "_"):
    link_sym = link_sym[:-4] + link_sym[-3:]
    if (len(link_sym) < 4): return None
  nnn = link_sym[:-3]
  try: nnn = int(nnn)
  except ValueError: return None
  if (nnn < 1 or nnn > len(symmetry_operators)): return None
  shifts = []
  for m in link_sym[-3:]:
    try: s = int(m)
    except ValueError: return None
    assert 0 <= s <= 9
    shifts.append(s-5)
  return symmetry_operators[nnn-1] + shifts
