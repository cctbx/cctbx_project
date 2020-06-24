
"""
Utility script to calculate per-residue RSCCs for a model versus an EM map with
an arbitrary origin.

Author: Nat Echols
Reference:
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS.
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Manuscript in review.
"""

from __future__ import absolute_import, division, print_function
import iotbx.phil
from cctbx import maptbx
from scitbx.array_family import flex
import sys

master_phil_str = """
model = None
  .type = path
map = None
  .type = path
d_min = 3.0
  .type = float
  .help = Optional cutoff resolution for computing F(calc). This will not \
    affect the dimensions of the ultimate FC map.
atom_radius = 1.5
  .type = float
"""

def run(args, out=sys.stdout):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    map_file_def="map",
    usage_string="""\
em_rscc.py model.pdb map.ccp4

%s""" % __doc__)
  params = cmdline.work.extract()
  assert (not None in [params.model, params.map])
  pdb_in = cmdline.get_file(params.model).file_object
  m = cmdline.get_file(params.map).file_object
  print("Input electron density map:", file=out)
  print("m.all()   :", m.map_data().all(), file=out)
  print("m.focus() :", m.map_data().focus(), file=out)
  print("m.origin():", m.map_data().origin(), file=out)
  print("m.nd()    :", m.map_data().nd(), file=out)
  print("m.size()  :", m.map_data().size(), file=out)
  print("m.focus_size_1d():", m.map_data().focus_size_1d(), file=out)
  print("m.is_0_based()   :", m.map_data().is_0_based(), file=out)
  print("map: min/max/mean:", flex.min(m.map_data()), flex.max(m.map_data()), flex.mean(m.map_data()), file=out)
  print("unit cell:", m.unit_cell().parameters(), file=out)
  symm = m.unit_cell_crystal_symmetry()
  xrs = pdb_in.input.xray_structure_simple(crystal_symmetry=symm)
  print("Setting up electron scattering table (d_min=%g)" % params.d_min, file=out)
  xrs.scattering_type_registry(
    d_min=params.d_min,
    table="electron")
  fc = xrs.structure_factors(d_min=params.d_min).f_calc()
  cg = maptbx.crystal_gridding(
    unit_cell=symm.unit_cell(),
    space_group_info=symm.space_group_info(),
    pre_determined_n_real=m.map_data().all())
  fc_map = fc.fft_map(
    crystal_gridding=cg).apply_sigma_scaling().real_map_unpadded()
  assert (fc_map.all() == fc_map.focus() == m.map_data().all())
  em_data = m.map_data()
  unit_cell_for_interpolation = m.grid_unit_cell()
  frac_matrix = unit_cell_for_interpolation.fractionalization_matrix()
  sites_cart = xrs.sites_cart()
  sites_frac = xrs.sites_frac()
  print("PER-RESIDUE CORRELATION:", file=out)
  for chain in pdb_in.hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      i_seqs = residue_group.atoms().extract_i_seq()
      values_em = flex.double()
      values_fc = flex.double()
      for i_seq in i_seqs :
        rho_em = maptbx.non_crystallographic_eight_point_interpolation(
          map=em_data,
          gridding_matrix=frac_matrix,
          site_cart=sites_cart[i_seq])
        rho_fc = fc_map.eight_point_interpolation(sites_frac[i_seq])
        values_em.append(rho_em)
        values_fc.append(rho_fc)
      cc = flex.linear_correlation(x=values_em, y=values_fc).coefficient()
      print(residue_group.id_str(), cc, file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
