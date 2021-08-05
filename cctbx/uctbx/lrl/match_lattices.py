from __future__ import division
import cctbx_uctbx_lrl_ext as ext

def get_input():
  import sys
  if len(sys.argv) > 1:
    with open(sys.argv[1]) as f:
      result = [l.strip() for l in f.readlines()]
  else:
    inp = sys.stdin.readline()
    result = []
    while inp.lower() != 'end':
      result.append(inp.strip())
      inp = sys.stdin.readline()
  return result

def selling_reduce(reader, mat=None):
  lc = ext.LatticeConverter
  if mat is None: mat = ext.MatS6()
  lat, cell = reader.GetLattice(), reader.GetCell()
  result = lc.SellingReduceCell(lat, cell, mat)
  return result


inputs = get_input()
input_readers = []
for inp in inputs:
  reader = ext.LRL_ReadLatticeData()
  reader.CellReader(inp)
  input_readers.append(reader)

mat_reference = ext.MatS6()
selling_reduce(input_readers[0], mat=mat_reference)
mat_reference = ext.MatS6.Inverse(mat_reference)

vLat = [] # These are the input lattices, Selling reduced
for reader in input_readers:
  vLat.append(ext.S6(selling_reduce(reader)))

lm = ext.LRL_LatticeMatcher()
lm.SetReferenceLattice(vLat[0])

vs = []
for lat in vLat:
  match = lm.MatchReference(lat)
  vs.append(mat_reference*match)

DC_ref = ext.DC(vLat[0])
for i, (s6, lat) in enumerate(zip(vs, vLat)):
  #distDC = ext.DC.DistanceBetween(DC_ref, ext.DC(lat))
  print(i, ' ', ext.LRL_Cell_Degrees(s6))
