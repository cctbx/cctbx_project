# Generate SHELX LATT and SYMM cards for a given space group.

import traceback

from cctbx import sgtbx

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("sgsymbol", "P1"),
              ("convention", "")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = form[key[0]].value.strip()
    else:
      inp.__dict__[key[0]] = key[1]
  return inp

def LATT_SYMM(space_group):
  lines = []
  l = lines.append
  Z = space_group.conventional_centring_type_symbol()
  Z_dict = {
    "P": 1,
    "I": 2,
    "R": 3,
    "F": 4,
    "A": 5,
    "B": 6,
    "C": 7,
  }
  try:
    LATT_N = Z_dict[Z]
  except:
    raise RuntimeError, "Error: Lattice type not supported by SHELX."
  # N must be made negative if the structure is non-centrosymmetric.
  if (space_group.is_centric()):
    if (not space_group.is_origin_centric()):
      raise RuntimeError, "Error: " \
        + " SHELX manual: If the structure is centrosymmetric, the" \
        + " origin MUST lie on a center of symmetry."
  else:
    LATT_N = -LATT_N;
  l("LATT %d" % (LATT_N,))
  # The operator x,y,z is always assumed, so MUST NOT be input.
  for i in xrange(1, space_group.n_smx()):
    l("SYMM %s" % (space_group(i).as_xyz(0, 0, "XYZ", ","),))
  return lines

def run(cctbx_url, inp):
  print "Content-type: text/html"
  print
  print "<pre>"
  try:
    space_group_info = sgtbx.space_group_info(
      symbol=inp.sgsymbol,
      table_id=inp.convention)
    space_group_info.show_summary()
    print
    cards = LATT_SYMM(space_group_info.group())
    for card in cards:
      print card
  except RuntimeError, e:
    print e
  except AssertionError:
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]
  print "</pre>"
