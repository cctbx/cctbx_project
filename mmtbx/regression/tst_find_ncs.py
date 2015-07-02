from __future__ import division

from libtbx.utils import null_out
from cStringIO import StringIO
from mmtbx.ncs import ncs

def remove_blank(text):
  return text.replace("\n"," ").replace(".0 "," ").replace(". "," ").replace(" ","")

text="""

Summary of NCS information
Wed Jul  1 10:54:17 2015
source_info ncs.ncs_spec

new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    3.2161   -7.6562  -13.0956
CHAIN A
RMSD 0
MATCHING 7
  RESSEQ 92:98

new_operator

rota_matrix    0.4525   -0.3091   -0.8365
rota_matrix    0.7737    0.6025    0.1959
rota_matrix    0.4434   -0.7359    0.5118
tran_orth    -0.1318    4.4033   -9.9054

center_orth   -9.2303   -5.9527   -6.7960
CHAIN R
RMSD 0.0007
MATCHING 7
  RESSEQ 92:98


"""

expected_text_group_specification=text

expected_text_resolve="""
new_ncs_group
rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    3.2161   -7.6562  -13.0956

rota_matrix    0.4525   -0.3091   -0.8365
rota_matrix    0.7737    0.6025    0.1959
rota_matrix    0.4434   -0.7359    0.5118
tran_orth    -0.1318    4.4033   -9.9054

center_orth   -9.2303   -5.9527   -6.7960

"""

expected_text_refine="""
refinement.pdb_interpretation.ncs_group {
  reference = chain 'A' and (resseq 92:98 )
  selection = chain 'R' and (resseq 92:98 )
}
"""

def tst_01():
  print "Read ncs-spec file and write out text...",

  f=open('ncs.ncs_spec','w')
  print >>f,text
  f.close()


  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.read_ncs('ncs.ncs_spec',quiet=True)

  f=StringIO()
  ncs_object.format_all_for_group_specification(out=f)
  found_text=f.getvalue()
  expected_text=expected_text_group_specification

  if remove_blank(" ".join(found_text.split("source_info")[1:]) ) != \
     remove_blank(" ".join(expected_text.split("source_info")[1:])):
    print "Expected: \n%s \nFound: \n%s" %(expected_text_group_specification,found_text)
    raise AssertionError, "FAILED"

  f=StringIO()
  ncs_object.format_all_for_resolve(out=f)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text_resolve):
    print "Expected: \n%s \nFound: \n%s" %(expected_text_resolve,found_text)
    raise AssertionError, "FAILED"

  f=StringIO()
  ncs_object.format_all_for_phenix_refine(out=f)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text_refine):
    print "Expected: \n%s \nFound: \n%s" %(expected_text_refine,found_text)
    raise AssertionError, "FAILED"

  print "OK"

if __name__=="__main__":
  tst_01()
