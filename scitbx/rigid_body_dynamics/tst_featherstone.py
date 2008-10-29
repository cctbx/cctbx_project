try:
  import scitbx
except ImportError:
  scitbx = None

if (scitbx is not None):
  from scitbx.rigid_body_dynamics import featherstone as fs
  from scitbx import matrix
  from libtbx.test_utils import approx_equal
else:
  import featherstone as fs
  import matrix
  def approx_equal(a1, a2): return True
  print "libtbx.test_utils not available: approx_equal() disabled"
  def sum(l):
    result = 0
    for e in l: result += e
    return result

import sys

def exercise_basic():
  assert approx_equal(sum(fs.Xrotx(0.1)), 5.98001666111)
  assert approx_equal(sum(fs.Xroty(0.2)), 5.92026631136)
  assert approx_equal(sum(fs.Xrotz(0.3)), 5.8213459565)
  assert approx_equal(sum(fs.Xtrans((1,2,3))), 6)
  assert approx_equal(sum(fs.crm((1,2,3,4,5,6))), 0)
  assert approx_equal(sum(fs.crf((1,2,3,4,5,6))), 0)
  assert approx_equal(sum(fs.mcI(
    m=1.234,
    c=matrix.col((1,2,3)),
    I=matrix.sym(sym_mat3=(2,3,4,0.1,0.2,0.3)))), 21.306)

def write_matlab_script(label, precision, code):
  open("tst_featherstone_exercise_%s.m" % label, "w").write("""\
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

def exercise_floatbase():
  write_matlab_script(label="floatbase", precision=7, code="""\
floatbase(autoTree(2))
""")
  tree = fs.floatbase(model=fs.autoTree(nb=2))
  assert tree.NB == 7
  assert tree.pitch == [fs.Inf, fs.Inf, fs.Inf, 0, 0, 0, 0]
  assert tree.parent == [-1, 0, 1, 2, 3, 4, 5]
  assert len(tree.Xtree) == 7
  assert len(tree.I) == 7
  expected = extract_sqr_matrices("""\
       0.0000000   0.0000000  -1.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       1.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000  -1.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   1.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000

       0.0000000   0.0000000   1.0000000   0.0000000   0.0000000   0.0000000
       1.0000000   0.0000000  -0.0000000   0.0000000   0.0000000   0.0000000
      -0.0000000   1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000  -0.0000000
       0.0000000   0.0000000   0.0000000  -0.0000000   1.0000000   0.0000000

       1.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   1.0000000   0.0000000   0.0000000   0.0000000
       0.0000000  -1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000
       0.0000000   0.0000000   0.0000000   0.0000000  -1.0000000   0.0000000

       0.0000000   0.0000000  -1.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       1.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000  -1.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   1.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000

       0.0000000   0.0000000   1.0000000   0.0000000   0.0000000   0.0000000
       1.0000000   0.0000000  -0.0000000   0.0000000   0.0000000   0.0000000
      -0.0000000   1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000  -0.0000000
       0.0000000   0.0000000   0.0000000  -0.0000000   1.0000000   0.0000000

       1.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   1.0000000   0.0000000   0.0000000   0.0000000
       0.0000000  -1.0000000   0.0000000   0.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000
       0.0000000   0.0000000   0.0000000   0.0000000  -1.0000000   0.0000000

       1   0   0   0   0   0
       0   1   0   0   0   0
       0   0   1   0   0   0
       0   0   0   1   0   0
       0   0   1   0   1   0
       0  -1   0   0   0   1
""")
  assert approx_equal(tree.Xtree, expected)
  for m in tree.I[:5]:
    assert approx_equal(m, [0]*36)
  assert approx_equal(tree.I[5:], expected_autoTree_I)

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

def exercise_IDf_FDf():
  write_matlab_script(label="IDf_FDf", precision=12, code="""\
% modified version of test_fb1.m

% test the correctness of the floating-base forward and inverse dynamics
% functions by checking that IDf is indeed the inverse of FDf.

% step 1: create a floating-base kinematic tree, and adjust some of the
% pitches so that it contains revolute, prismatic and helical joints.

NRJ = 8;                                % number of real joints
tree = floatbase( autoTree( NRJ+1, 2, 1, 0.9 ) );
tree.pitch(7) = 0.1;                    % 1st real joint
tree.pitch(8) = inf;                    % 2nd real joint

% step 2: select random initial conditions.

rand("state", [0]);
q   = pi * (2*rand(NRJ,1) - 1);
qd  = 2*rand(NRJ,1) - 1;
qdd = 2*rand(NRJ,1) - 1;
X = Xrotz(2*pi*rand(1)) * Xroty(2*pi*rand(1)) * Xrotx(2*pi*rand(1)) * ...
    Xtrans(2*rand(3,1)-1);
v = 2*rand(6,1) - 1;
q
qd
qdd
X
v

for i = 1:tree.NB
  f_ext{i} = 2*rand(6,1) - 1;
end
f_ext

grav_accn = 2*rand(3,1) - 1;
grav_accn

% step 3: use IDf to calculate the force required to produce qdd, and the
% resulting floating-base acceleration.  Then use FDf to calculate the
% accelerations produced by this force.

[a1,tau1] = IDf( tree, X, v, q, qd, qdd );
[a2,qdd2] = FDf( tree, X, v, q, qd, tau1 );

% step 4: compare results.  In theory, we should have a1==a2 and
% qdd==qdd2.  However, rounding errors will make them slightly different.
% Expect rounding errors in the vicinity of 1e-13 to 1e-14 on this test.

a1
tau1
a2
qdd2

[a1,tau1] = IDf( tree, X, v, q, qd, qdd, f_ext );
[a2,qdd2] = FDf( tree, X, v, q, qd, tau1, f_ext );

a1
tau1
a2
qdd2

[a1,tau1] = IDf( tree, X, v, q, qd, qdd, f_ext, grav_accn );
[a2,qdd2] = FDf( tree, X, v, q, qd, tau1, f_ext, grav_accn );

a1
tau1
a2
qdd2
""")
  NRJ = 8
  tree = fs.floatbase(fs.autoTree(nb=NRJ+1, bf=2, skew=1, taper=0.9))
  tree.pitch[6] = 0.1
  tree.pitch[7] = fs.Inf
  q = [
     2.16406631697e+00,
     1.62077531448e+00,
    -4.99063476296e-01,
    -1.51477073237e+00,
     7.08411636458e-02,
    -5.97316430786e-01,
     1.78315912482e+00,
    -1.23582258961e+00]
  qd = [
    -4.68060916953e-02,
     1.66764078910e-01,
     8.16225770391e-01,
     9.37371163478e-03,
    -4.36324311201e-01,
     5.11608408314e-01,
     2.36737993351e-01,
    -4.98987317275e-01]
  qdd = [
     8.19492511936e-01,
     9.65570952075e-01,
     6.20434471993e-01,
     8.04331900879e-01,
    -3.79704861361e-01,
     4.59663496520e-01,
     7.97676575936e-01,
     3.67967863831e-01]
  X = matrix.sqr([
    -7.94098097238e-01, -3.93496242589e-01, -4.63215844970e-01,
     0.00000000000e+00,  0.00000000000e+00,  0.00000000000e+00,
    -1.40429985987e-01,  8.60296737199e-01, -4.90070344951e-01,
     0.00000000000e+00,  0.00000000000e+00,  0.00000000000e+00,
     5.91343919389e-01, -3.24114533820e-01, -7.38418673903e-01,
     0.00000000000e+00,  0.00000000000e+00,  0.00000000000e+00,
    -1.54108230435e-02, -6.38333251412e-01,  5.68675368225e-01,
    -7.94098097238e-01, -3.93496242589e-01, -4.63215844970e-01,
    -1.20764881016e+00, -2.23662167230e-02,  3.06789675741e-01,
    -1.40429985987e-01,  8.60296737199e-01, -4.90070344951e-01,
    -3.07482337743e-01,  7.15611700425e-01, -5.60343309354e-01,
     5.91343919389e-01, -3.24114533820e-01, -7.38418673903e-01])
  v = matrix.col([
    -4.59804468946e-02,
     7.30619855543e-01,
    -4.79015379216e-01,
     6.10055654026e-01,
     9.73986076712e-02,
    -9.71916599672e-01])
  f_ext_given = [matrix.col(f) for f in [
    [],
    [],
    [],
    [],
    [],
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
     -3.34927701561e-01],
    [ 6.31826193067e-01,
     -7.98784959568e-01,
     -7.07283022175e-01,
      3.95341280382e-01,
     -9.09531864269e-01,
      1.47732073578e-01],
    [ 8.20032029398e-01,
      6.83959365214e-02,
      3.61178265125e-01,
     -9.46606410676e-01,
      2.69999819823e-01,
      2.12676835508e-01]]]
  grav_accn_given = [
     1.51905896063e-01,
    -2.17581181354e-01,
    -2.59720119330e-01]
  expected_a1_and_tau1 = [
    ([ 2.35908383275e-01,
       4.86085511739e-01,
      -1.92602724651e-01,
       1.80396138585e-01,
       7.70288216964e-01,
      -9.47954096322e+00],
     [ 3.42561006395e-01,
       1.67175023782e-01,
      -1.31504036908e-03,
       3.01430820237e-02,
      -1.70249903025e-02,
       4.07188292033e-04,
       2.50381946244e-02,
      -5.81335924462e-03]),
    ([ 2.59514051877e-01,
       1.66098102380e+00,
      -1.53687749224e+00,
       2.85174032290e-01,
      -1.47894556707e-01,
      -1.01591878867e+01],
     [ 3.09938773925e+00,
      -9.43607828571e-01,
      -5.13463254427e-01,
       2.08424562558e-01,
      -3.34937422780e-01,
       7.48377265553e-01,
       7.20070968322e-01,
      -3.47243185983e-01]),
    ([ 2.59514051877e-01,
       1.66098102380e+00,
      -1.53687749224e+00,
       4.37079928353e-01,
      -3.65475738061e-01,
      -6.08908005995e-01],
     [ 3.09938773925e+00,
      -9.43607828571e-01,
      -5.13463254427e-01,
       2.08424562558e-01,
      -3.34937422780e-01,
       7.48377265553e-01,
       7.20070968322e-01,
      -3.47243185983e-01])]
  for f_ext,grav_accn,(expected_a1,expected_tau1) in zip(
       [None, f_ext_given, f_ext_given],
       [None, None, grav_accn_given],
       expected_a1_and_tau1):
    a1, tau1 = fs.IDf(
      model=tree,
      Xfb=X, vfb=v,
      q=q, qd=qd, qdd=qdd,
      f_ext=f_ext, grav_accn=grav_accn)
    assert len(a1) == len(expected_a1)
    assert len(tau1) == len(expected_tau1)
    assert approx_equal(a1.elems, expected_a1)
    assert approx_equal([m.elems[0] for m in tau1], expected_tau1)
    continue
    a2, qdd2 = fs.FDf(
      model=tree,
      Xfb=X, vfb=v,
      q=q, qd=qd, tau=tau1,
      f_ext=f_ext, grav_accn=grav_accn)
    assert len(a2) == len(expected_a1)
    assert len(qdd2) == len(qdd)
    assert approx_equal(a2.elems, expected_a1)
    assert approx_equal([m.elems[0] for m in qdd2], qdd)

def exercise_standalone(tmpdir="tst_featherstone_tmpdir"):
  from libtbx.utils import copy_file
  import libtbx.load_env
  import os
  if (not os.path.isdir(tmpdir)):
    os.mkdir(tmpdir)
  os.chdir(tmpdir)
  scitbx_dist = libtbx.env.dist_path(module_name="scitbx")
  def cp(file_name):
    copy_file(os.path.join(scitbx_dist, file_name), ".")
  cp("matrix.py")
  cp("rigid_body_dynamics/featherstone.py")
  cp("rigid_body_dynamics/tst_featherstone.py")

def run(args):
  assert len(args) == 0
  exercise_basic()
  exercise_autoTree()
  exercise_floatbase()
  exercise_ID_FDab()
  exercise_IDf_FDf()
  if (scitbx is not None):
    exercise_standalone()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
