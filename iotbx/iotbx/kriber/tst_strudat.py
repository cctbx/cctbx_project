from iotbx.kriber import strudat
from cctbx import uctbx
from cStringIO import StringIO

def exercise():
  test_file = StringIO("""
*tric
Title
Reference
P 1
 11 12 13 100 110 120
Si 0.1 0.2 0.3
---------------------------------
*mono_b
Title
Reference
P 1 2 1
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*mono_c
Title
Reference
P 1 1 2
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*mono_a
Title
Reference
P 2 1 1
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*orth
Title
Reference
P 2 2 2
 11 12 13
Si 0.1 0.2 0.3 # remark
---------------------------------
*tetr
Title
Reference
P 4
 11 13
Si 0.1 0.2 0.3 4
---------------------------------
*trig
Title
Reference
R 3
 11 13
Si 0.1 0.2 0.3
---------------------------------
*rhom
Title
Reference
R 3 R
 11 100
Si 0.1 0.2 0.3
---------------------------------
*hexa
Title
Reference
P 6
 11 13
Si 0.1 0.2 0.3
---------------------------------
*cubi
Title
Reference
P 2 3
 11
Si 0.1 0.2 0.3
O  0.0 0.0 0.0
---------------------------------
""")
  all_entries = strudat.read_all_entries(test_file)
  for tag,cell in (("tric", (11,12,13,100,110,120)),
                   ("mono_b", (11,12,13,90,100,90)),
                   ("mono_c", (11,12,13,90,90,100)),
                   ("mono_a", (11,12,13,100,90,90)),
                   ("orth", (11,12,13,90,90,90)),
                   ("tetr", (11,11,13,90,90,90)),
                   ("trig", (11,11,13,90,90,120)),
                   ("rhom", (11,11,11,100,100,100)),
                   ("hexa", (11,11,13,90,90,120)),
                   ("cubi", (11,11,11,90,90,90))):
    assert all_entries.get(tag).unit_cell().is_similar_to(
      uctbx.unit_cell(cell))
  assert all_entries.get("orth").atoms[0].connectivity is None
  assert all_entries.get("tetr").atoms[0].connectivity == 4
  assert all_entries.get("cubi").as_xray_structure().scatterers().size() == 2
  print "OK"

if (__name__ == "__main__"):
  exercise()
