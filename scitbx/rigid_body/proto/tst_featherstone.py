from __future__ import absolute_import, division, print_function
from six.moves import zip
try:
  import scitbx
except ImportError:
  scitbx = None

if (scitbx is not None):
  from scitbx.rigid_body.proto import featherstone as fs
  from scitbx import matrix
  from libtbx.test_utils import approx_equal
else:
  from scitbx.rigid_body.proto import featherstone as fs
  import scitbx_matrix as matrix
  def approx_equal(a1, a2): return True
  print("libtbx.test_utils not available: approx_equal() disabled")
  def sum(l):
    result = 0
    for e in l: result += e
    return result

import sys

def exercise_basic():
  assert approx_equal(sum(fs.Xrotx(0.1)), 5.98001666111)
  assert approx_equal(sum(fs.Xroty(0.2)), 5.92026631136)
  assert approx_equal(sum(fs.Xrotz(0.3)), 5.8213459565)
  assert approx_equal(sum(fs.Xrot((1,2,3,4,5,6,7,8,9))), 90)
  assert approx_equal(sum(fs.Xtrans((1,2,3))), 6)
  assert approx_equal(sum(fs.crm((1,2,3,4,5,6))), 0)
  assert approx_equal(sum(fs.crf((1,2,3,4,5,6))), 0)
  assert approx_equal(sum(fs.mcI(
    m=1.234,
    c=matrix.col((1,2,3)),
    I=matrix.sym(sym_mat3=(2,3,4,0.1,0.2,0.3)))), 21.306)

def write_matlab_script(label, precision, code):
  with open("tst_featherstone_exercise_%s.m" % label, "w") as f:
    f.write("""\
split_long_rows(0);
output_precision(%d);

%s""" % (precision, code))

def extract_sqr_matrices(s):
  result = []
  elems = []
  for line in s.splitlines():
    if (len(s) == 0): continue
    elems.extend([float(v) for v in line.split()])
    if (len(elems) == 36):
      result.append(matrix.sqr(elems=elems))
      elems = []
  assert len(elems) == 0
  return result

expected_autoTree_I = extract_sqr_matrices("""
       0.0025000   0.0000000   0.0000000   0.0000000  -0.0000000   0.0000000
       0.0000000   0.3345833   0.0000000   0.0000000   0.0000000  -0.5000000
       0.0000000   0.0000000   0.3345833  -0.0000000   0.5000000   0.0000000
       0.0000000   0.0000000  -0.0000000   1.0000000   0.0000000   0.0000000
      -0.0000000   0.0000000   0.5000000   0.0000000   1.0000000   0.0000000
       0.0000000  -0.5000000   0.0000000   0.0000000   0.0000000   1.0000000

       0.0025000   0.0000000   0.0000000   0.0000000  -0.0000000   0.0000000
       0.0000000   0.3345833   0.0000000   0.0000000   0.0000000  -0.5000000
       0.0000000   0.0000000   0.3345833  -0.0000000   0.5000000   0.0000000
       0.0000000   0.0000000  -0.0000000   1.0000000   0.0000000   0.0000000
      -0.0000000   0.0000000   0.5000000   0.0000000   1.0000000   0.0000000
       0.0000000  -0.5000000   0.0000000   0.0000000   0.0000000   1.0000000
""")

def exercise_autoTree():
  write_matlab_script(label="autoTree", precision=7, code="""\
autoTree(2)
""")
  tree = fs.autoTree(nb=2)
  assert tree.NB == 2
  assert tree.pitch == [0, 0]
  assert tree.parent == [-1, 0]
  assert len(tree.Xtree) == 2
  assert len(tree.I) == 2
  expected = extract_sqr_matrices("""
       1   0   0   0   0   0
       0   1   0   0   0   0
       0   0   1   0   0   0
       0   0  -0   1   0   0
      -0   0   0   0   1   0
       0  -0   0   0   0   1

       1   0   0   0   0   0
       0   1   0   0   0   0
       0   0   1   0   0   0
       0   0   0   1   0   0
       0   0   1   0   1   0
       0  -1   0   0   0   1
""")
  assert approx_equal(tree.Xtree, expected)
  assert approx_equal(tree.I, expected_autoTree_I)

def exercise_ID_FDab():
  write_matlab_script(label="ID_FDab", precision=12, code="""\
% modified version of test_FD.m

% test the correctness of FDab, and ID by checking that FDab
% is the inverse of ID.

% step 1: create a complicated kinematic tree, and adjust some of the
% pitches so that it contains helical and prismatic as well as revolute
% joints.

tree = autoTree( 12, 1.5, 1, 0.95 );
tree.pitch(3) = 0.1;
tree.pitch(5) = inf;
tree.pitch(7) = -0.1;
tree.pitch(9) = inf;

% step 2: choose random initial conditions

rand("state", [0]);
q   = pi * (2*rand(12,1) - 1);
qd  = 2*rand(12,1) - 1;
qdd = 2*rand(12,1) - 1;
q
qd
qdd

for i = 1:tree.NB
  f_ext{i} = 2*rand(6,1) - 1;
end
f_ext

grav_accn = 2*rand(3,1) - 1;
grav_accn

% step 3: use ID to calculate the force required to produce qdd; then use
% FDab to calculate the acceleration that this force
% produces.

tau = ID( tree, q, qd, qdd );
qdd_ab = FDab( tree, q, qd, tau );

% step 4: compare results.  In theory, we should have qdd_ab==qdd and
% qdd_crb==qdd.  However, rounding errors will make them slightly
% different.  Expect rounding errors in the vicinity of 1e-14 on this
% test.

tau
qdd_ab

tau = ID( tree, q, qd, qdd, f_ext );
qdd_ab = FDab( tree, q, qd, tau, f_ext );

tau
qdd_ab

tau = ID( tree, q, qd, qdd, f_ext, grav_accn );
qdd_ab = FDab( tree, q, qd, tau, f_ext, grav_accn );

tau
qdd_ab
""")
  tree = fs.autoTree(nb=12, bf=1.5, skew=1, taper=0.95)
  tree.pitch[2] = 0.1
  tree.pitch[4] = fs.Inf
  tree.pitch[6] = -0.1
  tree.pitch[8] = fs.Inf
  q = [
     2.16406631697e+00,
     1.62077531448e+00,
    -4.99063476296e-01,
    -1.51477073237e+00,
     7.08411636458e-02,
    -5.97316430786e-01,
     1.78315912482e+00,
    -1.23582258961e+00,
    -1.47045673813e-01,
     5.23904805187e-01,
     2.56424888393e+00,
     2.94483836087e-02]
  qd = [
    -4.36324311201e-01,
     5.11608408314e-01,
     2.36737993351e-01,
    -4.98987317275e-01,
     8.19492511936e-01,
     9.65570952075e-01,
     6.20434471993e-01,
     8.04331900879e-01,
    -3.79704861361e-01,
     4.59663496520e-01,
     7.97676575936e-01,
     3.67967863831e-01]
  qdd = [
    -5.57145690946e-02,
    -7.98597583863e-01,
    -1.31656329092e-01,
     2.21773946888e-01,
     8.26022106476e-01,
     9.33212735542e-01,
    -4.59804468946e-02,
     7.30619855543e-01,
    -4.79015379216e-01,
     6.10055654026e-01,
     9.73986076712e-02,
    -9.71916599672e-01]
  f_ext_given = [matrix.col(f) for f in [
    [ 4.39409372808e-01,
     -2.02352915551e-01,
      6.49689954296e-01,
      3.36306402464e-01,
     -9.97714361371e-01,
     -1.28442670694e-02],
    [ 7.35205550986e-01,
     -5.12178246226e-01,
     -3.49591274505e-01,
      7.40942464217e-01,
     -6.17865816995e-01,
      1.35021481241e-01],
    [-5.22768142770e-01,
      9.35080500580e-01,
      6.06358938560e-01,
     -1.04060857129e-01,
     -8.39108362895e-01,
     -3.59890790655e-01],
    [ 1.58812850411e-02,
      8.65667648454e-01,
     -7.81884308138e-01,
      1.02534492181e-01,
      4.13122819734e-01,
      9.48818226568e-02],
    [ 6.28933726583e-01,
      8.05672139406e-02,
      9.27677091948e-01,
      2.06371255923e-01,
      1.75234128351e-01,
     -1.10021947449e-01],
    [ 1.92573723166e-01,
     -2.30197708055e-01,
      1.51302028330e-01,
     -4.19340995194e-01,
     -6.21217342891e-01,
     -6.26540943489e-01],
    [ 2.25546359737e-01,
      3.13318777979e-01,
     -4.69380159812e-02,
     -8.20351277609e-01,
      5.15207843933e-01,
      7.53540741646e-01],
    [ 8.46762031893e-01,
      6.84920446280e-01,
      7.96346242716e-01,
      8.46164879640e-01,
      8.11998498961e-02,
     -2.17407899531e-01],
    [ 4.10566799709e-01,
     -4.48731757376e-01,
      6.23257417016e-01,
      6.98971930373e-01,
      7.90077934853e-01,
      1.79602367062e-01],
    [ 8.99529746464e-01,
      1.59390021491e-01,
     -9.88737867377e-02,
      3.20490757245e-01,
      9.92515678707e-01,
      8.33882435895e-01],
    [ 5.86650168260e-01,
     -8.35254023607e-01,
      2.25566210081e-01,
     -2.71115960617e-02,
      2.60294680823e-01,
      6.90155151343e-01],
    [-5.13928755876e-01,
      4.62978441582e-01,
     -7.65731413583e-01,
     -5.59078926264e-01,
      5.89165943421e-01,
     -3.34927701561e-01]]]
  grav_accn_given = [
     6.31826193067e-01,
    -7.98784959568e-01,
    -7.07283022175e-01]
  expected_taus = [
    [-1.07802735301e+01,
     -2.49508384454e+00,
      1.26947454308e+01,
      1.91173545407e+01,
     -2.15095801697e+00,
      4.51666158929e+00,
      4.74401858265e-01,
      2.18443796607e+00,
     -2.16077117905e+00,
      2.83314247697e-01,
     -1.05175849728e-01,
     -1.50102709018e-01],
    [-9.58228539302e+00,
      1.26086473525e+00,
      1.17897084957e+01,
      2.39352891888e+01,
     -2.05738994237e+00,
      3.95457944509e+00,
      1.53558904623e-01,
      1.01233113762e+00,
     -2.34037354611e+00,
      3.82188034435e-01,
     -3.30742059809e-01,
      6.15628704565e-01],
    [-1.50057630860e+01,
     -1.28108933617e+01,
     -3.95574970014e+00,
      6.23340162074e+00,
      3.55385144046e+00,
      1.48006620138e+00,
     -7.33767304124e-01,
     -1.54512467666e-01,
     -8.05660270955e-01,
      3.27795875394e-01,
     -2.60889565195e-01,
      8.55164362101e-01]]
  for f_ext,grav_accn,expected_tau in zip(
       [None, f_ext_given, f_ext_given],
       [None, None, grav_accn_given],
       expected_taus):
    tau = fs.ID(
      model=tree, q=q, qd=qd, qdd=qdd, f_ext=f_ext, grav_accn=grav_accn)
    assert len(tau) == len(expected_tau)
    assert approx_equal([m.elems[0] for m in tau], expected_tau)
    qdd_ab = fs.FDab(
      model=tree, q=q, qd=qd, tau=tau, f_ext=f_ext, grav_accn=grav_accn)
    assert len(qdd_ab) == len(qdd)
    assert approx_equal([m.elems[0] for m in qdd_ab], qdd)

def exercise_standalone(tmpdir="tst_featherstone_tmpdir"):
  from libtbx.utils import copy_file
  import libtbx.load_env
  import os
  if (not os.path.isdir(tmpdir)):
    os.mkdir(tmpdir)
  os.chdir(tmpdir)
  scitbx_dist = libtbx.env.dist_path(module_name="scitbx")
  def cp(file_name, target="."):
    copy_file(os.path.join(scitbx_dist, file_name), target)
  cp("matrix/__init__.py", "scitbx_matrix.py")
  cp("rigid_body/proto/featherstone.py")
  cp("rigid_body/proto/tst_featherstone.py")

def run(args):
  assert len(args) == 0
  exercise_basic()
  exercise_autoTree()
  exercise_ID_FDab()
  if (scitbx is not None):
    exercise_standalone()
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
