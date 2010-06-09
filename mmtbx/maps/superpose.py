
from __future__ import division
import mmtbx.maps.utils
import iotbx.pdb
from iotbx import crystal_symmetry_from_any
from cctbx import maptbx
from cctbx.array_family import flex
from scitbx.math import superpose
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry, null_out
import libtbx.load_env
from libtbx import easy_run
from libtbx import adopt_init_args
import sys, os

def mask_grid(xrs, buffer, map_data, n_real):
  # XXX move to C++
  frac_min, frac_max = xrs.unit_cell().box_frac_around_sites(
    sites_cart = xrs.sites_cart(), buffer = buffer-1.5)
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  new_map = flex.double(flex.grid(n_real),0)
  for i in range(gridding_first[0], gridding_last[0]):
    for j in xrange(gridding_first[1], gridding_last[1]):
      for k in xrange(gridding_first[2], gridding_last[2]):
        if(i> 0 and i<n_real[0] and
           j> 0 and j<n_real[1] and
           k> 0 and k<n_real[2]):
          new_map[(i,j,k)] = map_data[(i,j,k)]
  return new_map

# this assumes that coordinates have already been shifted to all-positive
def generate_p1_box (pdb_hierarchy, buffer=10.0) :
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  xyz_max = sites_cart.max()
  a = xyz_max[0] + buffer
  b = xyz_max[1] + buffer
  c = xyz_max[2] + buffer
  combined = "%.3f,%.3f,%.3f,90,90,90,P1" % (a, b, c)
  symm = crystal_symmetry_from_any.from_string(combined)
  return symm

class common_frame_of_reference (object) :
  def __init__ (self, all_sites_cart, lsq_fits, buffer=10.0, log=sys.stdout) :
    fitted_sites = []
    original_sites = []
    minima = flex.vec3_double()
    for sites_cart, lsq_fit in zip(all_sites_cart, lsq_fits) :
      fitted_sites.append(sites_cart.deep_copy())
      minima.append(sites_cart.min())
      if lsq_fit is None :
        original_sites.append(sites_cart.deep_copy())
      else :
        old_sites = lsq_fit.r.inverse().elems * (sites_cart - lsq_fit.t.elems)
        original_sites.append(old_sites)
    xyz_min = minima.min()
    dxyz = (buffer - xyz_min[0], buffer - xyz_min[1], buffer - xyz_min[2])
    self.shifted_sites = []
    self.transformation_matrices = []
    for i, sites_cart in enumerate(fitted_sites) :
      new_sites_cart = sites_cart + dxyz
      #print new_sites_cart.min()
      self.shifted_sites.append(new_sites_cart)
      lsq_fit = superpose.least_squares_fit(
        reference_sites=new_sites_cart,
        other_sites=original_sites[i])
      self.transformation_matrices.append(lsq_fit.rt())

  def apply_new_sites_to_hierarchies (self, pdb_hierarchies) :
    assert len(pdb_hierarchies) == len(self.shifted_sites)
    for i, pdb_hierarchy in enumerate(pdb_hierarchies) :
      atoms = pdb_hierarchy.atoms()
      atoms.set_xyz(self.shifted_sites[i])

  def apply_new_sites_to_xray_structure (self, xray_structures) :
    pass

  def get_inverse_matrix (self, i) :
    assert i < len(self.transformation_matrices)
    rt = self.transformation_matrices[i]
    return rt.inverse()

  def inverse_transform_hierarchy (self, i, pdb_hierarchy) :
    rt = self.get_inverse_matrix(i)
    atoms = pdb_hierarchy.atoms()
    sites_cart = atoms.extract_xyz()
    atoms.set_xyz(rt.r.elems * sites_cart + rt.t.elems)
    uij = flex.sym_mat3_double(sites_cart.size(),[-1,-1,-1,-1,-1,-1])
    atoms.set_uij(uij)

def transform_map_by_lsq_fit (fft_map,
                              unit_cell,
                              lsq_fit_obj,
                              pdb_hierarchy,
                              d_min,
                              buffer=10.0,
                              resolution_factor=0.25,
                              file_name=None,
                              log=sys.stdout,
                              mask_map_grid=False,
                              format="xplor") :
  assert (format in ["xplor", "ccp4"])
  fake_symm = generate_p1_box(pdb_hierarchy, buffer=10.0)
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=fake_symm)
  f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
  fake_map = f_calc.fft_map(resolution_factor=resolution_factor)
  map_data_superposed = maptbx.superpose_maps(
    unit_cell_1        = unit_cell,
    unit_cell_2        = fake_symm.unit_cell(),
    map_data_1         = fft_map.real_map_unpadded(),
    n_real_2           = fake_map.n_real(),
    rotation_matrix    = lsq_fit_obj.r.elems,
    translation_vector = lsq_fit_obj.t.elems)
  if mask_map_grid :
    map_data_superposed = mask_grid(
      xrs      = xray_structure,
      buffer   = buffer,
      map_data = map_data_superposed,
      n_real   = fake_map.n_real())
  if file_name is not None :
    if format == "xplor" :
      print >> log, "    saving XPLOR map to %s" % file_name
      mmtbx.maps.utils.write_xplor_map(
        sites_cart=xray_structure.sites_cart(),
        unit_cell=fake_symm.unit_cell(),
        map_data=map_data_superposed,
        n_real=fake_map.n_real(),
        file_name=file_name,
        buffer=buffer)
    else :
      print >> log, "    saving CCP4 map to %s" % file_name
      mmtbx.maps.utils.write_ccp4_map(
        sites_cart=xray_structure.sites_cart(),
        unit_cell=fake_symm.unit_cell(),
        map_data=map_data_superposed,
        n_real=fake_map.n_real(),
        file_name=file_name,
        buffer=buffer)
  return xray_structure, map_data_superposed

# XXX: wrapper for multiprocessing
class transform_maps (object) :
  def __init__ (self,
                map_coeffs,
                map_types,
                pdb_hierarchy,
                unit_cell,
                lsq_fit_obj,
                output_files,
                resolution_factor=0.25,
                auto_run=True,
                format="ccp4") :
    adopt_init_args(self, locals())
    if auto_run :
      self.run()

  def run (self) :
    for (array, map_type, file_name) in zip(self.map_coeffs, self.map_types,
        self.output_files) :
      if array is None :
        continue
      (d_max, d_min) = array.d_max_min()
      fft_map = array.fft_map(resolution_factor=self.resolution_factor)
      fft_map.apply_sigma_scaling()
      transform_map_by_lsq_fit(
        fft_map=fft_map,
        unit_cell=self.unit_cell,
        lsq_fit_obj=self.lsq_fit_obj,
        pdb_hierarchy=self.pdb_hierarchy,
        d_min=d_min,
        file_name=file_name,
        log=null_out(),
        format=self.format)
      del fft_map

def exercise () :
  pdb_file_name_1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/1l3r_no_ligand.pdb",
    test=os.path.isfile)
  pdb_in_1 = iotbx.pdb.input(file_name = pdb_file_name_1)
  xrs1 = pdb_in_1.xray_structure_simple()
  ofn = "1l3r_rt.pdb"
  cmd = " ".join([
    "phenix.pdbtools",
    "%s"%pdb_file_name_1,
    "rotate='90 10 20' translate='10 10 10'",
    "output.file_name=%s"%ofn,
    "--quiet"])
  easy_run.call(cmd)
  pdb_in_rt = iotbx.pdb.input(file_name = ofn)
  xrs_rt = pdb_in_rt.xray_structure_simple()
  fft_map_1 = xrs1.structure_factors(d_min=1.5).f_calc().fft_map(
    resolution_factor = 1./3)
  fft_map_1.apply_sigma_scaling()
  map_data_1 = fft_map_1.real_map_unpadded()
  mmtbx.maps.utils.write_xplor_map(sites_cart = xrs1.sites_cart(),
    unit_cell  = xrs1.unit_cell(),
    map_data   = map_data_1,
    n_real     = fft_map_1.n_real(),
    file_name  = "1l3r.xplor")
  lsq_fit_obj = superpose.least_squares_fit(
    reference_sites=xrs_rt.sites_cart(),
    other_sites=xrs1.sites_cart().deep_copy())
  f_o_r = common_frame_of_reference(
    all_sites_cart=[xrs1.sites_cart(),xrs_rt.sites_cart()],
    lsq_fits=[None, lsq_fit_obj])
  hierarchy_rt = pdb_in_rt.construct_hierarchy()
  lsq_fit_obj = f_o_r.transformation_matrices[1]
  hierarchy_rt.atoms().set_xyz(f_o_r.shifted_sites[1])
  open("1l3r_rt.pdb", "w").write(hierarchy_rt.as_pdb_string())
  xrs_rt, map_data_rt = transform_map_by_lsq_fit(
    fft_map=fft_map_1,
    unit_cell=xrs1.unit_cell(),
    lsq_fit_obj=lsq_fit_obj.inverse(),
    pdb_hierarchy=hierarchy_rt,
    d_min=1.5,
    file_name="1l3r_rt.xplor",
    log=null_out())
  f_o_r.inverse_transform_hierarchy(1, hierarchy_rt)
  open("1l3r.pdb", "w").write(hierarchy_rt.as_pdb_string())
  #for sf1, sf2 in zip(xrs1.sites_frac(), xrs_rt.sites_frac()):
  #  e1 = map_data_1.eight_point_interpolation(sf1)
  #  e2 = map_data_rt.eight_point_interpolation(sf2)
  #  print abs(e1-e2)
  #  assert abs(e1-e2) < 1.
  print "OK"

if __name__ == "__main__" :
  exercise()
