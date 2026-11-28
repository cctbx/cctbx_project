from __future__ import absolute_import, division, print_function
import time
from mmtbx.refinement.real_space import fit_residues
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import crystal, uctbx
from scitbx.array_family import flex

"""
fit_residues.neighbors: atoms within a distance of the atoms of a residue.
It is called twice for every residue that is fitted. Done with the pair
generator over all atoms (unit cell and symmetry accounted for) it took 22 %
of side-chain fitting of a 113,000-atom model; if no periodic or symmetry
image can be within the distance (P1, the atoms in question far enough from
the faces of the cell: model in a box cut out of a cryo-EM map) plain
distances give the same atoms much faster.
"""

def get_sites(buffer_layer):
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = buffer_layer)
  ph.atoms().set_xyz(box.sites_cart)
  ph.atoms().reset_i_seq()
  return ph, box.crystal_symmetry()

def check(ph, crystal_symmetry, expect_fast):
  sps = crystal.special_position_settings(crystal_symmetry = crystal_symmetry)
  sites_cart = ph.atoms().extract_xyz()
  n_fast, n_all = 0, 0
  for radius in [1.5, 5, 10]:
    for rg in ph.residue_groups():
      i_selection = rg.atoms().extract_i_seq()
      general = fit_residues.neighbors_general(
        special_position_settings = sps, sites_cart = sites_cart,
        i_selection = i_selection, radius = radius)
      fast = fit_residues.neighbors_no_images(
        special_position_settings = sps, sites_cart = sites_cart,
        i_selection = i_selection, radius = radius)
      result = fit_residues.neighbors(
        special_position_settings = sps, sites_cart = sites_cart,
        i_selection = i_selection, radius = radius)
      assert list(result) == list(general)
      assert result.size() > i_selection.size()
      n_all += 1
      if(fast is not None):
        assert list(fast) == list(general)
        n_fast += 1
  print("fast way used: %d of %d"%(n_fast, n_all))
  if(expect_fast == "all"):    assert n_fast == n_all
  elif(expect_fast == "none"): assert n_fast == 0
  else:                        assert 0 < n_fast < n_all
  return n_fast

def exercise():
  # box much larger than the model: no images within reach
  ph, cs = get_sites(buffer_layer = 12)
  check(ph, cs, expect_fast = "all")
  # box that is just the model: images are within reach of atoms near faces,
  # they are found (the general way is used there)
  ph, cs = get_sites(buffer_layer = 3)
  check(ph, cs, expect_fast = "some")
  ph, cs = get_sites(buffer_layer = 0.5)
  check(ph, cs, expect_fast = "some")
  # symmetry: always the general way
  ph, cs = get_sites(buffer_layer = 12)
  cs = crystal.symmetry(unit_cell = cs.unit_cell(), space_group_symbol = "P 2 2 2")
  check(ph, cs, expect_fast = "none")
  # atoms outside of the unit cell: the general way
  ph, cs = get_sites(buffer_layer = 12)
  ph.atoms().set_xyz(ph.atoms().extract_xyz() + (200,0,0))
  check(ph, cs, expect_fast = "none")

def exercise_time():
  """
  Many atoms: the fast way is much faster.
  """
  ph, cs = get_sites(buffer_layer = 12)
  sps = crystal.special_position_settings(crystal_symmetry = cs)
  # 30000 atoms on a grid in a large box
  sites = flex.vec3_double()
  for i in range(30):
    for j in range(30):
      for k in range(34): sites.append((20+i*2., 20+j*2., 20+k*2.))
  cs = crystal.symmetry((100,100,110,90,90,90), "P1")
  sps = crystal.special_position_settings(crystal_symmetry = cs)
  i_selection = flex.size_t(range(15000, 15010))
  t0 = time.time()
  for it in range(3):
    r1 = fit_residues.neighbors_general(special_position_settings = sps,
      sites_cart = sites, i_selection = i_selection, radius = 7)
  t1 = time.time() - t0
  t0 = time.time()
  for it in range(3):
    r2 = fit_residues.neighbors(special_position_settings = sps,
      sites_cart = sites, i_selection = i_selection, radius = 7)
  t2 = time.time() - t0
  assert list(r1) == list(r2) and r1.size() > 100
  print("general %6.4f fast %6.4f s"%(t1, t2))
  assert t2 < t1/3, [t1, t2]

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  exercise_time()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
