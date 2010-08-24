from mmtbx.monomer_library import rna_sugar_pucker_analysis
from libtbx.test_utils import approx_equal
import sys

class atom(object):
  def __init__(O, xyz):
    O.xyz = xyz

def atom_dict(xyz_list):
  assert len(xyz_list) == 5
  return dict([(key, atom(xyz))
    for key,xyz in zip(["C1'", "O3'", "C3'", "C4'", "C5'"], xyz_list)])

def exercise(args):
  assert len(args) == 0
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  for i_pass,some_x in enumerate([31.101, 21.101]):
    analysis = rna_sugar_pucker_analysis.evaluate(
      params=params,
      residue_1_deoxy_ribo_atom_dict=atom_dict([
        (some_x,30.048,33.325), (22.184,28.704,36.539), (21.527,28.782,35.281),
        (22.379,28.287,34.119), (22.676,26.795,34.032)]),
      residue_1_c1p_outbound_atom=atom((19.628,30.157,32.953)),
      residue_2_p_atom=atom((33.058,33.903,32.613)))
    assert analysis.epsilon is None
    if (i_pass == 0):
      assert analysis.delta is None
    else:
      assert approx_equal(analysis.delta, 69.6819176834)
    assert analysis.p_distance_c1p_outbound_line is None
    assert analysis.is_epsilon_outlier is None
    if (i_pass == 0):
      assert analysis.is_2p_delta is None
      assert analysis.is_2p_p_distance_c1p_outbound_line is None
      assert analysis.is_2p_o3p_distance_c1p_outbound_line is None
      assert analysis.is_2p is None
    else:
      assert not analysis.is_2p_delta
      assert analysis.is_2p_p_distance_c1p_outbound_line is None
      assert not analysis.is_2p_o3p_distance_c1p_outbound_line
      assert not analysis.is_2p
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    residue_1_deoxy_ribo_atom_dict=atom_dict([
      (50.473,9.942,36.204), (48.952,10.236,33.644), (49.474,11.213,34.497),
      (50.978,11.174,34.278), (51.710,12.501,34.095)]),
    residue_1_c1p_outbound_atom=atom((50.926,10.084,37.569)),
    residue_2_p_atom=atom((48.114,10.641,32.347)))
  assert approx_equal(analysis.epsilon, 250.715932662)
  assert approx_equal(analysis.delta, 133.811229901)
  assert approx_equal(analysis.p_distance_c1p_outbound_line, 1.52374220064)
  assert approx_equal(analysis.o3p_distance_c1p_outbound_line, 0.860583842822)
  assert not analysis.is_epsilon_outlier
  assert analysis.is_2p_delta
  assert analysis.is_2p_p_distance_c1p_outbound_line
  assert analysis.is_2p_o3p_distance_c1p_outbound_line
  assert analysis.is_2p
  for residue_2_p_atom in [atom((41.856,17.581,31.191)), None]:
    analysis = rna_sugar_pucker_analysis.evaluate(
      params=params,
      residue_1_deoxy_ribo_atom_dict=atom_dict([
        (42.736,12.998,31.661), (42.922,16.529,31.706), (42.428,15.271,31.984),
        (43.132,14.710,33.185), (42.545,15.146,34.511)]),
      residue_1_c1p_outbound_atom=atom((41.433,12.332,31.466)),
      residue_2_p_atom=residue_2_p_atom)
    assert not analysis.is_epsilon_outlier
    assert not analysis.is_2p_delta
    if (residue_2_p_atom is not None):
      assert not analysis.is_2p_p_distance_c1p_outbound_line
    else:
      assert analysis.is_2p_p_distance_c1p_outbound_line is None
    assert not analysis.is_2p_o3p_distance_c1p_outbound_line
    assert not analysis.is_2p
  #
  params.epsilon_range_min = None
  params.epsilon_range_max = None
  params.delta_range_2p_min = None
  params.delta_range_2p_max = None
  params.delta_range_3p_min = None
  params.delta_range_3p_max = None
  params.p_distance_c1p_outbound_line_2p_max = None
  params.o3p_distance_c1p_outbound_line_2p_max = None
  analysis = rna_sugar_pucker_analysis.evaluate(
    params=params,
    residue_1_deoxy_ribo_atom_dict=atom_dict([
      (50.473,9.942,36.204), (48.952,10.236,33.644), (49.474,11.213,34.497),
      (50.978,11.174,34.278), (51.710,12.501,34.095)]),
    residue_1_c1p_outbound_atom=atom((50.926,10.084,37.569)),
    residue_2_p_atom=atom((48.114,10.641,32.347)))
  assert approx_equal(analysis.epsilon, 250.715932662)
  assert approx_equal(analysis.delta, 133.811229901)
  assert approx_equal(analysis.p_distance_c1p_outbound_line, 1.52374220064)
  assert analysis.is_epsilon_outlier is None
  assert analysis.is_2p_delta is None
  assert analysis.is_2p_p_distance_c1p_outbound_line is None
  assert analysis.is_2p_o3p_distance_c1p_outbound_line is None
  assert analysis.is_2p is None
  #
  print "OK"

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
