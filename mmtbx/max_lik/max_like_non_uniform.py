import math
from cctbx.array_family import flex
from cctbx.eltbx import van_der_waals_radii
from cctbx import maptbx
from cctbx import crystal
from libtbx.test_utils import approx_equal
from cctbx import miller
from scitbx import fftpack
from scitbx import matrix
from cctbx.eltbx.xray_scattering import wk1995

import boost.python
ext = boost.python.import_ext("mmtbx_max_lik_ext")
from mmtbx_max_lik_ext import *

import boost.python
ext = boost.python.import_ext("mmtbx_masks_ext")
from mmtbx_masks_ext import *



class ordered_solvent_distribution(object):

  def __init__(self,
               structure = None,
               fo = None,
               grid_step = None,
               rad=None,
               n_real = None,
               nshells = -1):
    assert [grid_step,n_real].count(None) == 1
    assert structure is not None
    assert fo is not None
    self.n_real = n_real
    self.nshells = nshells
    self._structure = structure
    self._grid_step = grid_step
    self._fo = fo
    self._f = None
    self._fc= None
    self._distribution = None
    if(rad is None):
       self.rad = 1.0
    else:
       self.rad = rad
    self._f_ordered_solvent()


  def _f_ordered_solvent(self):
    fo = self._fo
    fc = fo.structure_factors_from_scatterers(xray_structure = self._structure
                                             ).f_calc()
    self._fc = fc
    if(self.n_real is None):
       crystal_gridding = maptbx.crystal_gridding(unit_cell = self._structure.unit_cell(),
                                               step = self._grid_step)
       n_real = crystal_gridding.n_real()
    else:
       n_real = self.n_real

    xyzf = flex.vec3_double()
    atmrad = flex.double()
    elements = []
    for scatterer in self._structure.scatterers():
      xyzf.append(list(scatterer.site))
      atmrad.append(van_der_waals_radii.vdw.table[scatterer.element_symbol()])
      elements.append( scatterer.element_symbol() )
    assert xyzf.size() == atmrad.size()

# get residue name; wrong way to do the things ; must be improved later

    #sel_flag = flex.int(xyzf.size(),1)
    #if (getattr(self._structure, "scatterer_pdb_records", None) is not None):
    #  i=0
    #  for rec in self._structure.scatterer_pdb_records:
    #    if(rec.resName == "HOH"):
    #      sel_flag[i] = 0
    #    i+=1

    # !!! inverse logic for artificial case
    sel_flag = flex.int(xyzf.size(),1)
    if (getattr(self._structure, "scatterer_pdb_records", None) is not None):
      i=0
      for rec in self._structure.scatterer_pdb_records:
        if(rec.resName == "HOH"):
          sel_flag[i] = 0
        i+=1

    assert sel_flag.size() == atmrad.size()
    self._distribution = wat_dist()

    self._distribution.do_wat_dist(
          shell    = 0.0,
          xyzf     = xyzf,
          atmrad   = atmrad,
          element_symbol = elements,
          uc       = self._structure.unit_cell(),
          sg       = self._structure.space_group(),
          nxnynz   = n_real,
          sel_flag = sel_flag,
          rad      = self.rad,
          nshells  = self.nshells)
    data = self._distribution.data()

    ###############################
    #mask_data = mask(0.0,
    #                 1.0,
    #                 1.0,
    #                 xyzf,
    #                 atmrad,
    #                 self._structure.unit_cell(),
    #                 self._structure.space_group(),
    #                 crystal_gridding.n_real()).data()
    #
    #data.set_selected(mask_data == 0.0, 0.0)
    ###############################

    map_of_coeff = map_of_coeff_scaled(data,
                                       self._structure,
                                       n_real)

    from_map = maptbx.structure_factors.from_map(
       space_group=self._structure.space_group_info().group(),
       anomalous_flag=False,
       miller_indices=fo.indices(),
       complex_map=map_of_coeff,
       conjugate_flag=True)
    self._f = miller.array(
                miller_set = miller.set(
                   crystal_symmetry = crystal.symmetry(
                      unit_cell = self._structure.unit_cell(),
                      space_group_info = self._structure.space_group_info()),
                   indices = fo.indices(),
                   anomalous_flag = False),
                data = from_map.data())

    assert fc.indices().all_eq(self._f.indices()) == 1
    assert fc.indices().all_eq(fo.indices()) == 1
    assert flex.max(abs(self._f).data()) < 1.0
    return self

  def max_number_of_shells(self):
    return self._distribution.max_number_of_shells()

  def fcalc_from_distribution(self):
    return self._f

  #def f_ordered_solvent(self, n_water_atoms_absent = None,
  #                            bf_atoms_absent = None,
  #                            absent_atom_type = None):
  #  assert n_water_atoms_absent is not None
  #  nsym = self._f.space_group().order_z()
  #  n_lost_w = nsym * n_water_atoms_absent
  #  data = self._f.data() * n_lost_w
  #  d_star_sq_data = self._f.d_star_sq().data()
  #  table = wk1995(absent_atom_type).fetch()
  #  ff = table.at_d_star_sq(d_star_sq_data)
  #  factor = ff * flex.exp(-bf_atoms_absent/4.0*d_star_sq_data)
  #  f_by_m = miller.array(miller_set = self._f, data = data*factor)
  #  return f_by_m

  def distribution_as_array(self):
    return self._distribution.data()

  def distribution_as_xplor_file(self, file_name = None):
    if file_name is not None:
      self._distribution.as_xplor_map(self._structure.unit_cell(), file_name)

  def print_all_reflections_in_file(self, file_name = None):
    if(file_name is not None):
      file = open(file_name,"w")
      ic = self._fc.indices()
      i  = self._f .indices()
      io = self._fo.indices()
      fcd = abs(self._fc).data()
      fcp = self._fc.phases().data()
      fd  = abs(self._f).data()
      fp  = self._f.phases().data()
      fod = abs(self._fo).data()
      st = "%4d%4d%4d %10.6f %10.6f %4d%4d%4d %10.6f %10.6f %4d%4d%4d %10.3f \n"
      for j in xrange(fd.size()):
        file.write(st % (ic[j][0],ic[j][1],ic[j][2],fcd[j],fcp[j],\
                         i[j][0] ,i[j][1] ,i[j][2] ,fd[j] ,fp[j],\
                         io[j][0],io[j][1],io[j][2],fod[j]))

  #def fcalc_plus_f(self, n_water_atoms_absent = None):
  #  nsym = self._f.space_group().order_z()
  #  n_lost_w = nsym * n_water_atoms_absent
  #  #ss = 1./flex.pow2(self._f.d_spacings().data())
  #  data = self._fc.data() + self._f.data() * n_lost_w
  #  #data = self._fc.data().deep_copy()
  #  #for i in xrange(len(ss)):
  #  #  data[i] = self._fc.data()[i] + self._f.data()[i] * n_lost_w * form_factor("O", ss[i]) * \
  #  #          math.exp(-25.0/4.0*ss[i])
  #  fc_plus_f = miller.array(miller_set = self._fc,
  #                           data = data)
  #  return fc_plus_f

def form_factor(absent_atom_type,ss):
    #
    # W & K form-factor for an atom, B=0.0
    #
    table=wk1995(absent_atom_type).fetch()
    a_wk=table.array_of_a()
    b_wk=table.array_of_b()
    c_wk=table.c()
    result_wk=c_wk
    for i in xrange(5):
       result_wk += a_wk[i]*math.exp(-b_wk[i]*ss/4.0)
    return result_wk

def map_of_coeff_scaled(mask_map,structure,nxyz):
     assert mask_map.is_0_based()
     assert not mask_map.is_padded()
     fft_manager = fftpack.real_to_complex_3d(mask_map.focus())
     padded_data = maptbx.copy(mask_map,
                               flex.grid(fft_manager.m_real()
                              ).set_focus(fft_manager.n_real()))
     map_of_coeff = fft_manager.forward(padded_data)
     scale = matrix.col(nxyz).product()/structure.unit_cell().volume()
     sc = matrix.col(mask_map.focus()).product()/structure.unit_cell().volume()
     assert sc == scale
     map_of_coeff /= scale
     return map_of_coeff

def f_ordered_solvent(f,
                      n_water_atoms_absent,
                      bf_atoms_absent,
                      absent_atom_type):
  nsym = f.space_group().order_z()
  n_lost_w = nsym * n_water_atoms_absent
  data = f.data() * n_lost_w
  d_star_sq_data = f.d_star_sq().data()
  table = wk1995(absent_atom_type).fetch()
  ff = table.at_d_star_sq(d_star_sq_data)
  factor = ff * flex.exp(-bf_atoms_absent/4.0*d_star_sq_data)
  f_by_m = miller.array(miller_set = f, data = data*factor)
  return f_by_m

def alpha_beta(f_dist,
               n_atoms_included,
               n_nonwater_atoms_absent,
               n_water_atoms_absent,
               bf_atoms_absent,
               final_error,
               absent_atom_type):
  nsym = f_dist.space_group().order_z()
  ss = 1./flex.pow2(f_dist.d_spacings().data())
  n_part   = nsym * n_atoms_included
  n_lost_p = nsym * n_nonwater_atoms_absent
  n_lost_w = nsym * n_water_atoms_absent
  f_dist_data = flex.abs(f_dist.data())
  a_d = flex.exp( -0.25 * ss * final_error**2 * math.pi**3 )
  d_star_sq_data = f_dist.d_star_sq().data()
  assert approx_equal(ss,d_star_sq_data)
  table = wk1995(absent_atom_type).fetch()
  ff = table.at_d_star_sq(d_star_sq_data)
  factor = ff * flex.exp(-bf_atoms_absent/4.0*d_star_sq_data)
  b_d = ((1.-a_d*a_d)*n_part+n_lost_p+n_lost_w*(1.-f_dist_data*f_dist_data))*\
                                                                  factor*factor
  alpha = f_dist.array(data = a_d)
  beta  = f_dist.array(data = b_d)
  return alpha, beta
