from cctbx import crystal, sgtbx, xray
from cctbx import geometry_restraints
from cctbx.array_family import flex
import iotbx.cif.restraints
from libtbx.test_utils import show_diff

from cStringIO import StringIO

def exercise_geometry_restraints_as_cif():
  quartz = xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))
  bond_proxies = geometry_restraints.shared_bond_simple_proxy((
    geometry_restraints.bond_simple_proxy(
      i_seqs=[0,1],
      rt_mx_ji=sgtbx.rt_mx("x-y,x,z-2/3"),
      distance_ideal=1.6,
      weight=3.2),
    geometry_restraints.bond_simple_proxy(
      i_seqs=[0,1],
      distance_ideal=1.7,
      weight=1.8),
  ))
  angle_proxies = geometry_restraints.shared_angle_proxy((
    geometry_restraints.angle_proxy(
      i_seqs = [1,0,1],
      sym_ops = (sgtbx.rt_mx("x-y,x,z-2/3"),
                 sgtbx.rt_mx(),
                 sgtbx.rt_mx("-y,x-y,z-1/3")),
      angle_ideal=103,
      weight=2),
    geometry_restraints.angle_proxy(
      i_seqs = [1,0,1],
      sym_ops = (sgtbx.rt_mx("y+1,-x+y+1,z-1/3"),
                 sgtbx.rt_mx(),
                 sgtbx.rt_mx("-y,x-y,z-1/3")),
      angle_ideal=110,
      weight=5),
    geometry_restraints.angle_proxy(
      i_seqs = [0,1,0],
      sym_ops = (sgtbx.rt_mx("y,-x+y,z+2/3"),
                 sgtbx.rt_mx(),
                 sgtbx.rt_mx("-x+y,-x,z+1/3")),
      angle_ideal=150,
      weight=5),
  ))
  bond_similarity_proxies = geometry_restraints.shared_bond_similarity_proxy((
    geometry_restraints.bond_similarity_proxy(
      i_seqs=[(0,1),(0,1),(0,1)],
      sym_ops=(sgtbx.rt_mx("x-y,x,z-2/3"),
               sgtbx.rt_mx("-y,x-y,z-1/3"),
               sgtbx.rt_mx("y+1,-x+y+1,z-1/3")),
      weights=(1,1,1)),
  ))
  cif_block = iotbx.cif.model.block()
  iotbx.cif.restraints.add_to_cif_block(
    cif_block, quartz,
    bond_proxies=bond_proxies,
    angle_proxies=angle_proxies,
    bond_similarity_proxies=bond_similarity_proxies)
  s = StringIO()
  cif_block.show(out=s)
  assert not show_diff(s.getvalue(), """\
loop_
  _restr_distance_atom_site_label_1
  _restr_distance_atom_site_label_2
  _restr_distance_site_symmetry_2
  _restr_distance_target
  _restr_distance_target_weight_param
  _restr_distance_diff
   Si O 2_554 1.6000 0.5590 -0.0160
   Si O 1 1.7000 0.7454 -2.3838

loop_
  _restr_angle_atom_site_label_1
  _restr_angle_atom_site_label_2
  _restr_angle_atom_site_label_3
  _restr_angle_site_symmetry_1
  _restr_angle_site_symmetry_2
  _restr_angle_site_symmetry_3
  _restr_angle_target
  _restr_angle_target_weight_param
  _restr_angle_diff
   O Si O 1 2_554 1 103.0000 0.7071 1.6926
   O Si O 1 3_664 1 110.0000 0.4472 -1.3127
   Si O Si 3 1 3 150.0000 0.4472 3.0700

loop_
  _restr_equal_distance_class_class_id
  _restr_equal_distance_class_target_weight_param
  _restr_equal_distance_class_average
  _restr_equal_distance_class_esd
  _restr_equal_distance_class_diff_max
   1 1.0000 1.6160 0.0000 0.0000

loop_
  _restr_equal_distance_atom_site_label_1
  _restr_equal_distance_atom_site_label_2
  _restr_equal_distance_site_symmetry_2
  _restr_equal_distance_class_id
   Si O 2_554 1
   Si O 4_554 1
   Si O 3_664 1

""")

def exercise_adp_restraints_as_cif():
  import libtbx.load_env
  if not libtbx.env.has_module("smtbx"):
    print "Skipping exercise_adp_restraints_as_cif(): smtbx not available"
    return
  from smtbx.refinement.restraints import adp_restraints as smtbx_adp_restraints
  import smtbx.development
  xs = smtbx.development.sucrose()
  rigid_bond_proxies = smtbx_adp_restraints.rigid_bond_restraints(
    xray_structure=xs).proxies[:3]
  adp_similarity_proxies = smtbx_adp_restraints.adp_similarity_restraints(
    xray_structure=xs).proxies[:3]
  isotropic_adp_proxies = smtbx_adp_restraints.isotropic_adp_restraints(
    xray_structure=xs).proxies[:3]
  cif_block = iotbx.cif.model.block()
  iotbx.cif.restraints.add_to_cif_block(
    cif_block, xs,
    rigid_bond_proxies=rigid_bond_proxies,
    adp_similarity_proxies=adp_similarity_proxies,
    isotropic_adp_proxies=isotropic_adp_proxies)
  s = StringIO()
  cif_block.show(out=s)
  assert not show_diff(s.getvalue(), """\
loop_
  _restr_U_rigid_atom_site_label_1
  _restr_U_rigid_atom_site_label_2
  _restr_U_rigid_target_weight_param
  _restr_U_rigid_U_parallel
  _restr_U_rigid_diff
   O1 C1 0.0100 0.0176 0.0006
   O1 C2 0.0100 0.0194 -0.0053
   O1 C3 0.0100 0.0177 0.0013

loop_
  _restr_U_similar_atom_site_label_1
  _restr_U_similar_atom_site_label_2
  _restr_U_similar_weight_param
   O1 C1 0.0400
   O1 C6 0.0400
   O2 C2 0.0800

loop_
  _restr_U_iso_atom_site_label
  _restr_U_iso_weight_param
   O1 0.1000
   O2 0.2000
   O3 0.2000

""")


def run():
  exercise_adp_restraints_as_cif()
  exercise_geometry_restraints_as_cif()
  print "OK"

if __name__ == '__main__':
  run()
