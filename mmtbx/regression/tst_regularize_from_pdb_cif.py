from __future__ import absolute_import, division, print_function

from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
from mmtbx.secondary_structure.regularize_from_pdb import \
   replace_with_segments_from_pdb

def remove_blank(text):
  return text.replace(" ","").replace("\n","")

pdb_str_1="""
ATOM      1  CA  GLY U 111       6.802  13.378  -2.270  1.00 52.86           C
ATOM      2  CA  GLY U 112       4.244  15.964  -0.444  1.00 52.86           C
ATOM      3  CA  GLY U 113       1.213  15.704  -1.022  1.00 52.86           C
ATOM      4  CA  GLY U 114      -0.519  14.243  -3.109  1.00 52.86           C
ATOM      5  CA  GLY U 115      -2.874  13.497  -4.987  1.00 52.86           C
ATOM      6  CA  GLY U 116      -0.713  13.888  -7.199  1.00 52.86           C
ATOM      7  CA  GLY U 117       1.212  11.376  -6.966  1.00 52.86           C
ATOM      8  CA  GLY U 118      -1.819  10.012  -8.645  1.00 52.86           C
ATOM      9  CA  GLY U 119       0.375   9.344 -10.927  1.00 52.86           C
ATOM     10  CA  GLY U 120       2.317   5.692 -10.865  1.00 52.86           C
ATOM     11  CA  GLY U 121       4.732   3.967 -11.502  1.00 52.86           C
ATOM     12  CA  GLY U 122       8.166   5.698 -10.793  1.00 42.55           C
ATOM     13  CA  GLY U 123       6.434   6.101  -6.390  1.00 52.86           C
ATOM     14  CA  GLY U 124       9.700   6.022  -5.333  1.00 52.86           C
ATOM     15  CA  GLY U 125       8.739   8.281  -3.579  1.00 52.86           C
ATOM     16  CA  GLY U 126       6.548   8.326  -1.298  1.00 52.86           C
ATOM     17  CA  GLY U 127       4.110   9.670  -0.100  1.00 52.86           C
ATOM     18  CA  GLY U 128       3.406  11.974   2.546  1.00 52.86           C
ATOM     19  CA  GLY U 129       0.747  12.768   3.800  1.00 52.86           C
ATOM     20  CA  GLY U 130      -1.068  13.851   6.552  1.00 52.86           C
ATOM     21  CA  GLY U 131      -3.058  16.387   5.476  1.00 52.86           C
ATOM     22  CA  GLY U 132      -4.637  14.251   3.493  1.00 52.86           C
ATOM     23  CA  GLY U 133      -4.850  12.995   0.227  1.00 52.86           C
ATOM     24  CA  GLY U 134      -7.881  13.254  -1.451  1.00 52.86           C
ATOM     25  CA  GLY U 135     -10.306  12.656  -3.983  1.00 52.86           C
TER
END
"""

# Convert to mmcif:
chain_addition = "ZXLONG"
from libtbx.test_utils import convert_pdb_to_cif_for_pdb_str
convert_pdb_to_cif_for_pdb_str(locals(),chain_addition=chain_addition)

def tst_01():
  print("Regularizing allowing insertions...", end=' ')
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(pdb_str_1)).construct_hierarchy()
  r=replace_with_segments_from_pdb(args=[],pdb_hierarchy=hierarchy,
    out=null_out())

  expected_text="""
ID: 1 ChainID: 'U%s'  RMSD:  1.31 A  (n=22) Junction RMSD:  0.68 A (n=7)
Complete: True  Insertions/deletions: True
Input model start: 111  end: 135  length: 25
Replacement start: 111  end: 136  length: 26
""" %(chain_addition)
  f=StringIO()
  for rss in r.model_replacement_segment_summaries:
    rss.show_summary(out=f)

  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("Expected: \n%s \nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")
  print("OK")

def tst_02():
  print("Regularizing not allowing insertions...", end=' ')
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(pdb_str_1)).construct_hierarchy()
  r=replace_with_segments_from_pdb(args=['alpha.allow_insertions=false'],
   pdb_hierarchy=hierarchy,
    out=null_out())

  expected_text="""
ID: 1 ChainID: 'U%s'  RMSD:  1.28 A  (n=25) Junction RMSD:  0.68 A (n=8)
Complete: True  Insertions/deletions: False
Input model start: 111  end: 135  length: 25
Replacement start: 111  end: 135  length: 25
""" %(chain_addition)
  f=StringIO()
  for rss in r.model_replacement_segment_summaries:
    rss.show_summary(out=f)

  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("Expected: \n%s \nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")
  print("OK")


if __name__=="__main__":
  tst_01()
  tst_02()
