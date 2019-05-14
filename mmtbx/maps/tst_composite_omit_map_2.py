
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import iotbx.pdb
import mmtbx.f_model
from cctbx import miller
import mmtbx.maps
from cctbx import maptbx
import mmtbx.maps.composite_omit_map
import sys

pdb_str="""\n
REMARK iotbx.pdb.box_around_molecule --buffer-layer=2 "m.pdb"
REMARK Date 2013-12-16 Time 13:47:17 PST -0800 (1387230437.01 s)
CRYST1   12.380   14.615   11.502  90.00  90.00  90.00 P 1
ATOM      1  CB  PHE A   1       4.641   4.458   4.502  1.00 10.00           C
ATOM      2  CG  PHE A   1       3.809   3.637   5.453  1.00 10.00           C
ATOM      3  CD1 PHE A   1       2.792   2.781   4.971  1.00 10.00           C
ATOM      4  CD2 PHE A   1       4.035   3.712   6.843  1.00 10.00           C
ATOM      5  CE1 PHE A   1       2.000   2.000   5.869  1.00 10.00           C
ATOM      6  CE2 PHE A   1       3.256   2.941   7.761  1.00 10.00           C
ATOM      7  CZ  PHE A   1       2.234   2.081   7.270  1.00 10.00           C
ATOM      8  C   PHE A   1       4.830   6.416   2.964  1.00 10.00           C
ATOM      9  O   PHE A   1       5.380   5.842   2.000  1.00 10.00           O
ATOM     10  OXT PHE A   1       5.017   7.615   3.259  1.00 10.00           O
ATOM     11  N   PHE A   1       2.749   5.066   3.014  1.00 10.00           N
ATOM     12  CA  PHE A   1       3.874   5.605   3.831  1.00 10.00           C
ATOM     13  CB  PHE A   2       9.641   9.458   9.502  1.00 10.00           C
ATOM     14  C   PHE A   2       9.830  11.416   7.964  1.00 10.00           C
ATOM     15  O   PHE A   2      10.380  10.842   7.000  1.00 10.00           O
ATOM     16  OXT PHE A   2      10.017  12.615   8.259  1.00 10.00           O
ATOM     17  N   PHE A   2       7.749  10.066   8.014  1.00 10.00           N
ATOM     18  CA  PHE A   2       8.874  10.605   8.831  1.00 10.00           C
TER
END
"""

def get_cc(mc1, mc2, xrs):
  crystal_gridding = mc1.crystal_gridding(
    d_min             = mc1.d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 0.25)
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = mc1)
  fft_map.apply_sigma_scaling()
  m1 = fft_map.real_map_unpadded()
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = mc2)
  fft_map.apply_sigma_scaling()
  m2 = fft_map.real_map_unpadded()
  assert m1.focus()==m2.focus()
  assert m1.all()==m2.all()
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = mc1.unit_cell(),
    fft_n_real = m1.focus(),
    fft_m_real = m1.all(),
    sites_cart = flex.vec3_double(xrs.sites_cart()),
    site_radii = flex.double([1.5]*xrs.scatterers().size()))
  cc = flex.linear_correlation(x=m1.select(sel), y=m2.select(sel)).coefficient()
  def md(m, xrs):
    r = flex.double()
    for sf in xrs.sites_frac():
      r.append(m.eight_point_interpolation(sf))
    return flex.mean(r)
  return cc, md(m=m1, xrs=xrs), md(m=m2, xrs=xrs)

def run(args):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph = pdb_inp.construct_hierarchy()
  ph.write_pdb_file(file_name="m.pdb")
  xrs = pdb_inp.xray_structure_simple()
  sel = ph.atom_selection_cache().selection("resseq 1")
  xrs1 = xrs.select( sel)
  xrs2 = xrs.select(~sel)
  #
  F  =  xrs.structure_factors(d_min=1.5).f_calc()
  FB = xrs1.structure_factors(d_min=1.5).f_calc()
  FC =   FB.phase_transfer(phase_source=F)
  #
  mtz_dataset = abs(FB).as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(miller_array=FB.generate_r_free_flags(),
    column_root_label="R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "data.mtz")
  #
  mtz_dataset = F.as_mtz_dataset(column_root_label="f")
  mtz_dataset.add_miller_array(miller_array=FB, column_root_label="FB")
  mtz_dataset.add_miller_array(miller_array=FC, column_root_label="FC")
  mtz_dataset.add_miller_array(miller_array=abs(FB), column_root_label="FOBS")
  #
  fmodel = mmtbx.f_model.manager(f_obs = abs(FB), xray_structure = xrs)
  # mc is model biased: one should not see any density for resseq 2
  mc = fmodel.electron_density_map().map_coefficients(
    map_type     = "Fo",
    isotropize   = False,
    fill_missing = False)
  #
  print(get_cc(mc1=mc, mc2=F,  xrs=xrs1), get_cc(mc1=mc, mc2=F,  xrs=xrs2))
  print(get_cc(mc1=mc, mc2=FB, xrs=xrs1), get_cc(mc1=mc, mc2=FB, xrs=xrs2))
  print(get_cc(mc1=mc, mc2=FC, xrs=xrs1), get_cc(mc1=mc, mc2=FC, xrs=xrs2))
  print()
  #
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = fmodel.f_obs().d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 0.25)
  #
  oo = mmtbx.maps.composite_omit_map.run(
    crystal_gridding     = crystal_gridding,
    full_resolution_map  = False,
    fmodel               = fmodel.deep_copy(),
    box_size_as_fraction = 0.2)
  mco = oo.map_coefficients(filter_noise=False)
  print(get_cc(mc1=mco, mc2=F, xrs=xrs1), get_cc(mc1=mco, mc2=F, xrs=xrs2))
  #
  mtz_dataset.add_miller_array(miller_array=mco, column_root_label="O")
  mtz_dataset.add_miller_array(miller_array=mc, column_root_label="mc")
  #
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "m.mtz")

if (__name__ == "__main__"):
  run(sys.argv[1:])
