import iotbx.cns.xray_structure
import iotbx.cns.miller_array
import iotbx.cns.reflection_reader
import iotbx.cns.crystal_symmetry_utils
import iotbx.cns.crystal_symmetry_from_inp
from iotbx.cns import sdb_reader
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from cStringIO import StringIO
import re
import sys

def exercise_crystal_symmetry_utils():
  as_rem = iotbx.cns.crystal_symmetry_utils.crystal_symmetry_as_sg_uc
  as_inp = iotbx.cns.crystal_symmetry_utils.crystal_symmetry_as_cns_inp_defines
  crystal_symmetry = crystal.symmetry(unit_cell=None, space_group_info=None)
  sg_uc = as_rem(crystal_symmetry=crystal_symmetry)
  assert sg_uc == "sg=None a=None b=None c=None alpha=None beta=None gamma=None"
  sg_uc = as_inp(crystal_symmetry=crystal_symmetry)
  assert not show_diff("\n".join(sg_uc), """\
{===>} sg="None";
{===>} a=None;
{===>} b=None;
{===>} c=None;
{===>} alpha=None;
{===>} beta=None;
{===>} gamma=None;""")
  #
  crystal_symmetry = crystal.symmetry(
    unit_cell=(3,4,5,89,87,93),
    space_group_info=None)
  sg_uc = as_rem(crystal_symmetry=crystal_symmetry)
  assert sg_uc == "sg=None a=3 b=4 c=5 alpha=89 beta=87 gamma=93"
  sg_uc = as_inp(crystal_symmetry=crystal_symmetry)
  assert not show_diff("\n".join(sg_uc), """\
{===>} sg="None";
{===>} a=3;
{===>} b=4;
{===>} c=5;
{===>} alpha=89;
{===>} beta=87;
{===>} gamma=93;""")
  #
  crystal_symmetry = crystal.symmetry(
    unit_cell=None,
    space_group_symbol=19)
  sg_uc = as_rem(crystal_symmetry=crystal_symmetry)
  assert sg_uc == \
    "sg=P2(1)2(1)2(1) a=None b=None c=None alpha=None beta=None gamma=None"
  sg_uc = as_inp(crystal_symmetry=crystal_symmetry)
  assert not show_diff("\n".join(sg_uc), """\
{===>} sg="P2(1)2(1)2(1)";
{===>} a=None;
{===>} b=None;
{===>} c=None;
{===>} alpha=None;
{===>} beta=None;
{===>} gamma=None;""")
  #
  for symbols in sgtbx.space_group_symbol_iterator():
    cs1 = sgtbx.space_group_info(
      group=sgtbx.space_group(symbols)).any_compatible_crystal_symmetry(
        volume=1000)
    sg_uc = as_rem(crystal_symmetry=cs1)
    m = re.match(iotbx.cns.crystal_symmetry_utils.re_sg_uc, sg_uc)
    assert m is not None
    cs2 = iotbx.cns.crystal_symmetry_utils.crystal_symmetry_from_re_match(m=m)
    assert cs1.is_similar_symmetry(cs2)
    #
    sg_uc = as_inp(crystal_symmetry=cs1)
    sio = StringIO()
    print >> sio, "\n".join(sg_uc)
    sio = StringIO(sio.getvalue())
    cs2 = iotbx.cns.crystal_symmetry_from_inp.extract_from(file=sio)
    assert cs1.is_similar_symmetry(cs2)
  #
  uc_sg = """\
a= 82.901 b= 82.901 c= 364.175 alpha= 90 beta= 90 gamma= 120 sg= P6(5)22"""
  m = re.match(iotbx.cns.crystal_symmetry_utils.re_uc_sg, uc_sg)
  cs = iotbx.cns.crystal_symmetry_utils.crystal_symmetry_from_re_match(
    m=m, i_uc=1, i_sg=7)
  assert cs is not None
  assert str(cs.unit_cell()) == "(82.901, 82.901, 364.175, 90, 90, 120)"
  assert str(cs.space_group().info()) == "P 65 2 2"

def exercise_sdb(verbose=0):
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=False,
    random_u_iso=True)
  f_abs = abs(structure.structure_factors(
    anomalous_flag=False, d_min=2, algorithm="direct").f_calc())
  sdb_out = structure.as_cns_sdb_file(
    file="foo.sdb",
    description="random_structure",
    comment=["any", "thing"],
    group="best")
  if (0 or verbose):
    sys.stdout.write(sdb_out)
  sdb_files = sdb_reader.multi_sdb_parser(StringIO(sdb_out))
  assert len(sdb_files) == 1
  structure_read = sdb_files[0].as_xray_structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=structure.unit_cell(),
      space_group_info=None))
  f_read = abs(f_abs.structure_factors_from_scatterers(
    xray_structure=structure_read, algorithm="direct").f_calc())
  regression = flex.linear_regression(f_abs.data(), f_read.data())
  assert regression.is_well_defined()
  if (0 or verbose):
    regression.show_summary()
  assert abs(regression.slope()-1) < 1.e-4
  assert abs(regression.y_intercept()) < 1.e-3

def exercise_reflection_reader():
  crf = iotbx.cns.reflection_reader.cns_reflection_file
  try:
    # just to make sure a bug in handling {} doesn't get reintroduced
    crf(file_handle=StringIO("}{"))
  except iotbx.cns.reflection_reader.CNS_input_Error, e:
    assert str(e) == "premature end-of-file"
  else:
    raise Exception_expected
  #
  def check(expected):
    c = crf(file_handle=si)
    so = StringIO()
    c.crystal_symmetry().show_summary(f=so)
    assert not show_diff(so.getvalue(), expected)
  si = StringIO("""\
remark a= 40.000 b= 50.000 c=  60.000 alpha= 90 beta= 90 gamma= 90 sg= P2
 remark symop (X,Y,Z)
remark symop (-X,-Y,Z)
CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1 21 1
DECLare NAME=FOBS                   DOMAin=RECIprocal   TYPE=REAL END
INDE     1    2    3 FOBS=   380.500
""")
  check("""\
Unit cell: (40, 50, 60, 90, 90, 90)
Space group: P 1 1 2 (No. 3)
""")
  si = StringIO("""\
 remark a= 40.000 b= 50.000 c=  60.000 alpha= 90 beta= 90 gamma= 90 sg= P2
CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1 21 1
DECLare NAME=FOBS                   DOMAin=RECIprocal   TYPE=REAL END
INDE     1    2    3 FOBS=   380.500
""")
  check("""\
Unit cell: (40, 50, 60, 90, 90, 90)
Space group: P 1 2 1 (No. 3)
""")
  si = StringIO("""\
CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1 21 1
DECLare NAME=FOBS                   DOMAin=RECIprocal   TYPE=REAL END
INDE     1    2    3 FOBS=   380.500
""")
  check("""\
Unit cell: (10, 20, 30, 90, 90, 90)
Space group: P 1 21 1 (No. 4)
""")

def exercise_miller_array_as_cns_hkl():
  s = StringIO()
  crystal_symmetry = crystal.symmetry()
  for anomalous_flag in [False, True]:
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=flex.miller_index([(1,2,3),(-3,5,-7)]),
      anomalous_flag=anomalous_flag)
    for data in [flex.bool((False,True)),
                 flex.int((-3,4)),
                 flex.double((10,13)),
                 flex.complex_double((10,13)),
                 flex.hendrickson_lattman(((0,1,2,3),(4,5,6,7)))]:
      miller_array = miller_set.array(data=data)
      miller_array.export_as_cns_hkl(file_object=s)
      if (isinstance(data, flex.double)):
        miller_array = miller_set.array(data=data, sigmas=data/10.)
        miller_array.export_as_cns_hkl(file_object=s)
  assert not show_diff(s.getvalue(), """\
NREFlections=2
ANOMalous=FALSe
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=INTEger END
INDEx 1 2 3 DATA= 0
INDEx -3 5 -7 DATA= 1
NREFlections=2
ANOMalous=FALSe
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=INTEger END
INDEx 1 2 3 DATA= -3
INDEx -3 5 -7 DATA= 4
NREFlections=2
ANOMalous=FALSe
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=REAL END
INDEx 1 2 3 DATA= 10
INDEx -3 5 -7 DATA= 13
NREFlections=2
ANOMalous=FALSe
DECLare NAME=FOBS DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=SIGMA DOMAin=RECIprocal TYPE=REAL END
INDEx 1 2 3 FOBS= 10 SIGMA= 1
INDEx -3 5 -7 FOBS= 13 SIGMA= 1.3
NREFlections=2
ANOMalous=FALSe
DECLare NAME=F  DOMAin=RECIprocal TYPE=COMPLEX END
INDEx 1 2 3 F= 10 0
INDEx -3 5 -7 F= 13 0
NREFlections=2
ANOMalous=FALSe
DECLare NAME=PA  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PB  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PC  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PD  DOMAin=RECIprocal TYPE=REAL END
GROUp  TYPE=HL
     OBJEct=PA
     OBJEct=PB
     OBJEct=PC
     OBJEct=PD
END
INDEx 1 2 3 PA= 0 PB= 1 PC= 2 PD= 3
INDEx -3 5 -7 PA= 4 PB= 5 PC= 6 PD= 7
NREFlections=2
ANOMalous=TRUE
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=INTEger END
INDEx 1 2 3 DATA= 0
INDEx -3 5 -7 DATA= 1
NREFlections=2
ANOMalous=TRUE
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=INTEger END
INDEx 1 2 3 DATA= -3
INDEx -3 5 -7 DATA= 4
NREFlections=2
ANOMalous=TRUE
DECLare NAME=DATA  DOMAin=RECIprocal TYPE=REAL END
INDEx 1 2 3 DATA= 10
INDEx -3 5 -7 DATA= 13
NREFlections=2
ANOMalous=TRUE
DECLare NAME=FOBS DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=SIGMA DOMAin=RECIprocal TYPE=REAL END
INDEx 1 2 3 FOBS= 10 SIGMA= 1
INDEx -3 5 -7 FOBS= 13 SIGMA= 1.3
NREFlections=2
ANOMalous=TRUE
DECLare NAME=F  DOMAin=RECIprocal TYPE=COMPLEX END
INDEx 1 2 3 F= 10 0
INDEx -3 5 -7 F= 13 0
NREFlections=2
ANOMalous=TRUE
DECLare NAME=PA  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PB  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PC  DOMAin=RECIprocal TYPE=REAL END
DECLare NAME=PD  DOMAin=RECIprocal TYPE=REAL END
GROUp  TYPE=HL
     OBJEct=PA
     OBJEct=PB
     OBJEct=PC
     OBJEct=PD
END
INDEx 1 2 3 PA= 0 PB= 1 PC= 2 PD= 3
INDEx -3 5 -7 PA= 4 PB= 5 PC= 6 PD= 7
""")

def exercise_reflection_file_as_miller_array():
  refl_1 = """\
 NREFlection=         3
 ANOMalous=FALSe { equiv. to HERMitian=TRUE}
 DECLare NAME=FOBS         DOMAin=RECIprocal   TYPE=COMP END
 DECLare NAME=PHASE        DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=SIGMA        DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=TEST         DOMAin=RECIprocal   TYPE=INTE END
 DECLare NAME=FOM          DOMAin=RECIprocal   TYPE=REAL END
 INDE     1    1    2 FOBS=   713.650 PHASE=     0.000 SIGMA=     3.280
                   TEST=  0 FOM=    -1.000
 INDE    -2    0    3 FOBS=   539.520 PHASE=     0.000 SIGMA=     4.140
                   TEST=  1 FOM=    -1.000
 INDE     0    2    1 FOBS=   268.140 PHASE=     0.000 SIGMA=     1.690
                   TEST=  0 FOM=    -1.000
"""
  refl_2 = """\
 { sg=C2 a=97.37 b=46.64 c=65.47 alpha=90 beta=115.4 gamma=90 }
 NREFlection=         3
 ANOMalous=FALSe { equiv. to HERMitian=TRUE}
 DECLare NAME=FOBS         DOMAin=RECIprocal   TYPE=COMP END
 DECLare NAME=SIGMA        DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=TEST         DOMAin=RECIprocal   TYPE=INTE END
 DECLare NAME=FOM          DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=PA           DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=PB           DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=PC           DOMAin=RECIprocal   TYPE=REAL END
 DECLare NAME=PD           DOMAin=RECIprocal   TYPE=REAL END
 GROUp TYPE=HL
     OBJEct=PA
     OBJEct=PB
     OBJEct=PC
     OBJEct=PD
 END
 INDE     1    1    2 FOBS=   713.650    44.854 SIGMA=     3.280 TEST=         0
                   FOM=     0.044 PA=     0.066 PB=     0.065 PC=    -0.001
                   PD=    -0.099
 INDE    -2    0    3 FOBS=   539.520     0.000 SIGMA=     4.140 TEST=         1
                   FOM=     0.994 PA=     2.893 PB=     0.000 PC=     0.000
                   PD=     0.000
 INDE     0    2    1 FOBS=   268.140   184.247 SIGMA=     1.690 TEST=         0
                   FOM=     0.890 PA=   -46.660 PB=    -3.465 PC=   -12.988
                   PD=    -1.940
"""
  crystal_symmetry = crystal.symmetry(
    unit_cell=(97.37, 46.64, 65.47, 90, 115.4, 90),
    space_group_symbol="C 1 2 1")
  all_arrays = [iotbx.cns.reflection_reader.cns_reflection_file(
    file_handle=StringIO(refl)).as_miller_arrays(
      crystal_symmetry=cs)
        for refl,cs in [(refl_1,crystal_symmetry), (refl_2,None)]]
  for miller_arrays in all_arrays:
    for miller_array in miller_arrays:
      assert miller_array.crystal_symmetry().is_similar_symmetry(
        other=crystal_symmetry,
        relative_length_tolerance=1.e-10,
        absolute_angle_tolerance=1.e-10)
      assert not miller_array.anomalous_flag()
      assert list(miller_array.indices()) == [(1,1,2), (-2,0,3), (0,2,1)]
      lbl = miller_array.info().label_string()
      if (lbl == "FOBS,SIGMA"):
        assert str(miller_array.observation_type()) == "xray.amplitude"
        assert approx_equal(miller_array.data(), [713.65, 539.52, 268.14])
        assert approx_equal(miller_array.sigmas(), [3.28, 4.14, 1.69])
      else:
        assert miller_array.observation_type() is None
        if (lbl == "TEST"):
          assert list(miller_array.data()) == [0, 1, 0]
        elif (lbl == "PHASE"):
          assert approx_equal(miller_array.data(), [0, 0, 0])
        elif (lbl == "FOM"):
          if (miller_array.data()[0] < 0):
            assert approx_equal(miller_array.data(), [-1, -1, -1])
          else:
            assert approx_equal(miller_array.data(), [0.044, 0.994, 0.890])
        elif (lbl == "PA,PB,PC,PD"):
          assert approx_equal(miller_array.data(), [
            (0.066, 0.065, -0.001, -0.099),
            (2.893, 0.000, 0.000, 0.000),
            (-46.660, -3.465, -12.988, -1.940)])
        else:
          raise RuntimeError

def run():
  verbose = "--Verbose" in sys.argv[1:]
  exercise_crystal_symmetry_utils()
  exercise_sdb(verbose)
  exercise_reflection_reader()
  exercise_reflection_file_as_miller_array()
  exercise_miller_array_as_cns_hkl()
  print "OK"

if (__name__ == "__main__"):
  run()
