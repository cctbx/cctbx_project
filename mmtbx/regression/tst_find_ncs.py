from __future__ import division

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

text_helical="""
REMARK 350   BIOMT1   1 -0.982450 -0.186524  0.000000        0.00000
REMARK 350   BIOMT2   1  0.186524 -0.982450  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000      -35.25000
REMARK 350   BIOMT1   2 -0.980683  0.195604  0.000000        0.00000
REMARK 350   BIOMT2   2 -0.195604 -0.980683  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000      -33.84000
REMARK 350   BIOMT1   3 -0.835712  0.549169  0.000000        0.00000
REMARK 350   BIOMT2   3 -0.549169 -0.835712  0.000000        0.00000
REMARK 350   BIOMT3   3  0.000000  0.000000  1.000000      -32.43000
REMARK 350   BIOMT1   4 -0.568705  0.822541  0.000000        0.00000
REMARK 350   BIOMT2   4 -0.822541 -0.568705  0.000000        0.00000
REMARK 350   BIOMT3   4  0.000000  0.000000  1.000000      -31.02000
REMARK 350   BIOMT1   5 -0.218654  0.975802  0.000000        0.00000
REMARK 350   BIOMT2   5 -0.975802 -0.218654  0.000000        0.00000
REMARK 350   BIOMT3   5  0.000000  0.000000  1.000000      -29.61000
REMARK 350   BIOMT1   6  0.163326  0.986572  0.000000        0.00000
REMARK 350   BIOMT2   6 -0.986572  0.163326  0.000000        0.00000
REMARK 350   BIOMT3   6  0.000000  0.000000  1.000000      -28.20000
REMARK 350   BIOMT1   7  0.521456  0.853278  0.000000        0.00000
REMARK 350   BIOMT2   7 -0.853278  0.521456  0.000000        0.00000
REMARK 350   BIOMT3   7  0.000000  0.000000  1.000000      -26.79000
REMARK 350   BIOMT1   8  0.803441  0.595384  0.000000        0.00000
REMARK 350   BIOMT2   8 -0.595384  0.803441  0.000000        0.00000
REMARK 350   BIOMT3   8  0.000000  0.000000  1.000000      -25.38000
REMARK 350   BIOMT1   9  0.968104  0.250549  0.000000        0.00000
REMARK 350   BIOMT2   9 -0.250549  0.968104  0.000000        0.00000
REMARK 350   BIOMT3   9  0.000000  0.000000  1.000000      -23.97000
REMARK 350   BIOMT1  10  0.991399 -0.130872  0.000000        0.00000
REMARK 350   BIOMT2  10  0.130872  0.991399  0.000000        0.00000
REMARK 350   BIOMT3  10  0.000000  0.000000  1.000000      -22.56000
REMARK 350   BIOMT1  11  0.869926 -0.493183  0.000000        0.00000
REMARK 350   BIOMT2  11  0.493183  0.869926  0.000000        0.00000
REMARK 350   BIOMT3  11  0.000000  0.000000  1.000000      -21.15000
REMARK 350   BIOMT1  12  0.621421 -0.783477  0.000000        0.00000
REMARK 350   BIOMT2  12  0.783477  0.621421  0.000000        0.00000
REMARK 350   BIOMT3  12  0.000000  0.000000  1.000000      -19.74000
REMARK 350   BIOMT1  13  0.282174 -0.959363  0.000000        0.00000
REMARK 350   BIOMT2  13  0.959363  0.282174  0.000000        0.00000
REMARK 350   BIOMT3  13  0.000000  0.000000  1.000000      -18.33000
REMARK 350   BIOMT1  14 -0.098278 -0.995159  0.000000        0.00000
REMARK 350   BIOMT2  14  0.995159 -0.098278  0.000000        0.00000
REMARK 350   BIOMT3  14  0.000000  0.000000  1.000000      -16.92000
REMARK 350   BIOMT1  15 -0.464378 -0.885637  0.000000        0.00000
REMARK 350   BIOMT2  15  0.885637 -0.464378  0.000000        0.00000
REMARK 350   BIOMT3  15  0.000000  0.000000  1.000000      -15.51000
REMARK 350   BIOMT1  16 -0.762668 -0.646790  0.000000        0.00000
REMARK 350   BIOMT2  16  0.646790 -0.762668  0.000000        0.00000
REMARK 350   BIOMT3  16  0.000000  0.000000  1.000000      -14.10000
REMARK 350   BIOMT1  17 -0.949590 -0.313495  0.000000        0.00000
REMARK 350   BIOMT2  17  0.313495 -0.949590  0.000000        0.00000
REMARK 350   BIOMT3  17  0.000000  0.000000  1.000000      -12.69000
REMARK 350   BIOMT1  18 -0.997847  0.065577  0.000000        0.00000
REMARK 350   BIOMT2  18 -0.065577 -0.997847  0.000000        0.00000
REMARK 350   BIOMT3  18  0.000000  0.000000  1.000000      -11.28000
REMARK 350   BIOMT1  19 -0.900395  0.435074  0.000000        0.00000
REMARK 350   BIOMT2  19 -0.435074 -0.900395  0.000000        0.00000
REMARK 350   BIOMT3  19  0.000000  0.000000  1.000000       -9.87000
REMARK 350   BIOMT1  20 -0.671462  0.741039  0.000000        0.00000
REMARK 350   BIOMT2  20 -0.741039 -0.671462  0.000000        0.00000
REMARK 350   BIOMT3  20  0.000000  0.000000  1.000000       -8.46000
REMARK 350   BIOMT1  21 -0.344479  0.938794  0.000000        0.00000
REMARK 350   BIOMT2  21 -0.938794 -0.344479  0.000000        0.00000
REMARK 350   BIOMT3  21  0.000000  0.000000  1.000000       -7.05000
REMARK 350   BIOMT1  22  0.032806  0.999462  0.000000        0.00000
REMARK 350   BIOMT2  22 -0.999462  0.032806  0.000000        0.00000
REMARK 350   BIOMT3  22  0.000000  0.000000  1.000000       -5.64000
REMARK 350   BIOMT1  23  0.405301  0.914183  0.000000        0.00000
REMARK 350   BIOMT2  23 -0.914183  0.405301  0.000000        0.00000
REMARK 350   BIOMT3  23  0.000000  0.000000  1.000000       -4.23000
REMARK 350   BIOMT1  24  0.718612  0.695411  0.000000        0.00000
REMARK 350   BIOMT2  24 -0.695411  0.718612  0.000000        0.00000
REMARK 350   BIOMT3  24  0.000000  0.000000  1.000000       -2.82000
REMARK 350   BIOMT1  25  0.926988  0.375092  0.000000        0.00000
REMARK 350   BIOMT2  25 -0.375092  0.926988  0.000000        0.00000
REMARK 350   BIOMT3  25  0.000000  0.000000  1.000000       -1.41000
REMARK 350   BIOMT1  26  1.000000 -0.000000  0.000000        0.00000
REMARK 350   BIOMT2  26  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3  26  0.000000  0.000000  1.000000       -0.00000
REMARK 350   BIOMT1  27  0.926988 -0.375092  0.000000        0.00000
REMARK 350   BIOMT2  27  0.375092  0.926988  0.000000        0.00000
REMARK 350   BIOMT3  27  0.000000  0.000000  1.000000        1.41000
REMARK 350   BIOMT1  28  0.718612 -0.695411  0.000000        0.00000
REMARK 350   BIOMT2  28  0.695411  0.718612  0.000000        0.00000
REMARK 350   BIOMT3  28  0.000000  0.000000  1.000000        2.82000
REMARK 350   BIOMT1  29  0.405301 -0.914183  0.000000        0.00000
REMARK 350   BIOMT2  29  0.914183  0.405301  0.000000        0.00000
REMARK 350   BIOMT3  29  0.000000  0.000000  1.000000        4.23000
REMARK 350   BIOMT1  30  0.032806 -0.999462  0.000000        0.00000
REMARK 350   BIOMT2  30  0.999462  0.032806  0.000000        0.00000
REMARK 350   BIOMT3  30  0.000000  0.000000  1.000000        5.64000
REMARK 350   BIOMT1  31 -0.344479 -0.938794  0.000000        0.00000
REMARK 350   BIOMT2  31  0.938794 -0.344479  0.000000        0.00000
REMARK 350   BIOMT3  31  0.000000  0.000000  1.000000        7.05000
REMARK 350   BIOMT1  32 -0.671462 -0.741039  0.000000        0.00000
REMARK 350   BIOMT2  32  0.741039 -0.671462  0.000000        0.00000
REMARK 350   BIOMT3  32  0.000000  0.000000  1.000000        8.46000
REMARK 350   BIOMT1  33 -0.900395 -0.435074  0.000000        0.00000
REMARK 350   BIOMT2  33  0.435074 -0.900395  0.000000        0.00000
REMARK 350   BIOMT3  33  0.000000  0.000000  1.000000        9.87000
REMARK 350   BIOMT1  34 -0.997847 -0.065577  0.000000        0.00000
REMARK 350   BIOMT2  34  0.065577 -0.997847  0.000000        0.00000
REMARK 350   BIOMT3  34  0.000000  0.000000  1.000000       11.28000
REMARK 350   BIOMT1  35 -0.949590  0.313495  0.000000        0.00000
REMARK 350   BIOMT2  35 -0.313495 -0.949590  0.000000        0.00000
REMARK 350   BIOMT3  35  0.000000  0.000000  1.000000       12.69000
REMARK 350   BIOMT1  36 -0.762668  0.646790  0.000000        0.00000
REMARK 350   BIOMT2  36 -0.646790 -0.762668  0.000000        0.00000
REMARK 350   BIOMT3  36  0.000000  0.000000  1.000000       14.10000
REMARK 350   BIOMT1  37 -0.464378  0.885637  0.000000        0.00000
REMARK 350   BIOMT2  37 -0.885637 -0.464378  0.000000        0.00000
REMARK 350   BIOMT3  37  0.000000  0.000000  1.000000       15.51000
REMARK 350   BIOMT1  38 -0.098278  0.995159  0.000000        0.00000
REMARK 350   BIOMT2  38 -0.995159 -0.098278  0.000000        0.00000
REMARK 350   BIOMT3  38  0.000000  0.000000  1.000000       16.92000
REMARK 350   BIOMT1  39  0.282174  0.959363  0.000000        0.00000
REMARK 350   BIOMT2  39 -0.959363  0.282174  0.000000        0.00000
REMARK 350   BIOMT3  39  0.000000  0.000000  1.000000       18.33000
REMARK 350   BIOMT1  40  0.621421  0.783477  0.000000        0.00000
REMARK 350   BIOMT2  40 -0.783477  0.621421  0.000000        0.00000
REMARK 350   BIOMT3  40  0.000000  0.000000  1.000000       19.74000
REMARK 350   BIOMT1  41  0.869926  0.493183  0.000000        0.00000
REMARK 350   BIOMT2  41 -0.493183  0.869926  0.000000        0.00000
REMARK 350   BIOMT3  41  0.000000  0.000000  1.000000       21.15000
REMARK 350   BIOMT1  42  0.991399  0.130872  0.000000        0.00000
REMARK 350   BIOMT2  42 -0.130872  0.991399  0.000000        0.00000
REMARK 350   BIOMT3  42  0.000000  0.000000  1.000000       22.56000
REMARK 350   BIOMT1  43  0.968104 -0.250549  0.000000        0.00000
REMARK 350   BIOMT2  43  0.250549  0.968104  0.000000        0.00000
REMARK 350   BIOMT3  43  0.000000  0.000000  1.000000       23.97000
REMARK 350   BIOMT1  44  0.803441 -0.595384  0.000000        0.00000
REMARK 350   BIOMT2  44  0.595384  0.803441  0.000000        0.00000
REMARK 350   BIOMT3  44  0.000000  0.000000  1.000000       25.38000
REMARK 350   BIOMT1  45  0.521456 -0.853278  0.000000        0.00000
REMARK 350   BIOMT2  45  0.853278  0.521456  0.000000        0.00000
REMARK 350   BIOMT3  45  0.000000  0.000000  1.000000       26.79000
REMARK 350   BIOMT1  46  0.163326 -0.986572  0.000000        0.00000
REMARK 350   BIOMT2  46  0.986572  0.163326  0.000000        0.00000
REMARK 350   BIOMT3  46  0.000000  0.000000  1.000000       28.20000
REMARK 350   BIOMT1  47 -0.218654 -0.975802  0.000000        0.00000
REMARK 350   BIOMT2  47  0.975802 -0.218654  0.000000        0.00000
REMARK 350   BIOMT3  47  0.000000  0.000000  1.000000       29.61000
REMARK 350   BIOMT1  48 -0.568705 -0.822541  0.000000        0.00000
REMARK 350   BIOMT2  48  0.822541 -0.568705  0.000000        0.00000
REMARK 350   BIOMT3  48  0.000000  0.000000  1.000000       31.02000
REMARK 350   BIOMT1  49 -0.835712 -0.549169  0.000000        0.00000
REMARK 350   BIOMT2  49  0.549169 -0.835712  0.000000        0.00000
REMARK 350   BIOMT3  49  0.000000  0.000000  1.000000       32.43000
REMARK 350   BIOMT1  50 -0.980683 -0.195604  0.000000        0.00000
REMARK 350   BIOMT2  50  0.195604 -0.980683  0.000000        0.00000
REMARK 350   BIOMT3  50  0.000000  0.000000  1.000000       33.84000
"""
text_point_group="""
new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   18.1267   39.0869   16.3485

new_operator

rota_matrix    0.6235    0.7818    0.0000
rota_matrix   -0.7818    0.6235    0.0000
rota_matrix   -0.0000   -0.0000    1.0000
tran_orth   -14.8728   23.9186    0.0000

center_orth    8.7170   35.2568   16.3490

new_operator

rota_matrix   -0.2225    0.9749    0.0000
rota_matrix   -0.9749   -0.2225    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth    -5.4432   50.4629   -0.0007

center_orth    5.8467   25.5100   16.3495

new_operator

rota_matrix   -0.9010    0.4339   -0.0000
rota_matrix   -0.4339   -0.9010    0.0000
rota_matrix   -0.0000    0.0000    1.0000
tran_orth    21.1858   59.6326    0.0000

center_orth   11.6703   17.1839   16.3485

new_operator

rota_matrix   -0.9010   -0.4339   -0.0000
rota_matrix    0.4339   -0.9010   -0.0000
rota_matrix   -0.0000   -0.0000    1.0000
tran_orth    44.9594   44.5364   -0.0006

center_orth   21.8108   16.5521   16.3495

new_operator

rota_matrix   -0.2225   -0.9749    0.0000
rota_matrix    0.9749   -0.2225   -0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth    47.9895   16.5293    0.0005

center_orth   28.6368   24.0954   16.3473

new_operator

rota_matrix    0.6235   -0.7818    0.0000
rota_matrix    0.7818    0.6235    0.0000
rota_matrix   -0.0000    0.0000    1.0000
tran_orth    27.9737   -3.2911    0.0000

center_orth   26.9919   34.1217   16.3485
"""
text_nothing="""
new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   18.1267   39.0869   16.3485

new_operator

rota_matrix    0.6235    0.7818    0.0000
rota_matrix   -0.7818    0.6235    0.0000
rota_matrix   -0.0000   -0.0000    1.0000
tran_orth   -14.8728   23.9186    0.0000

center_orth    8.7170   35.2568   16.3490

new_operator

rota_matrix   -0.2225    0.9749    0.0000
rota_matrix   -0.9749   -0.2225    0.0000
rota_matrix    0.0000   -0.0000    1.0000
tran_orth    -5.4432   50.4629   -0.0007

center_orth    5.8467   25.5100   16.3495

new_operator

rota_matrix   -0.9010    0.4339   -0.0000
rota_matrix   -0.4339   -0.9010    0.0000
rota_matrix   -0.0000    0.0000    1.0000
tran_orth    21.1858   59.6326    0.0000

center_orth   11.6703   17.1839   16.3485

new_operator

rota_matrix   -0.9010   -0.4339   -0.0000
rota_matrix    0.4339   -0.9010   -0.0000
rota_matrix   -0.0000   -0.0000    1.0000
tran_orth    44.9594   44.5364   -0.0006

center_orth   21.8108   16.5521   16.3495

new_operator

rota_matrix   -0.2225   -0.9749    0.0000
rota_matrix    0.9749   -0.2225   -0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth    47.9895   16.5293    0.0005

center_orth   28.6368   24.0954   16.3473

"""
def tst_02():

  print "Test helical and point_group symmetry...",
  from mmtbx.ncs.ncs import ncs
  f=open('helical.ncs_spec','w')
  print >>f,text_helical
  f.close()
  ncs_object=ncs()
  ncs_object.read_ncs('helical.ncs_spec',quiet=True)
  assert ncs_object.is_helical_symmetry()
  assert not ncs_object.is_point_group_symmetry()

  f=open('point_group.ncs_spec','w')
  print >>f,text_point_group
  f.close()
  ncs_object=ncs()
  ncs_object.read_ncs('point_group.ncs_spec',quiet=True)
  assert not ncs_object.is_helical_symmetry()
  assert ncs_object.is_point_group_symmetry()

  f=open('nothing.ncs_spec','w')
  print >>f,text_nothing
  f.close()
  ncs_object=ncs()
  ncs_object.read_ncs('nothing.ncs_spec',quiet=True)
  assert not ncs_object.is_helical_symmetry()
  assert not ncs_object.is_point_group_symmetry()

  print "OK"

if __name__=="__main__":
  tst_01()
  tst_02()
