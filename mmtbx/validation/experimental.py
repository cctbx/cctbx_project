
# NOT USED ANYWHERE. STILL EXISTS FOR BACKWARDS COMPATIBILITY. PLEASE REMOVE ONCE READY!

"""
Model validation against experimental data, in both real and reciprocal space.
This does not actually handle any of the scaling and fmodel calculations,
which are performed approximately as in model_vs_data.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation import residue, validation
from libtbx import Auto, slots_getstate_setstate
from libtbx.str_utils import format_value
from libtbx.utils import null_out, Sorry
import sys

__real_space_attr__ = [
  "b_iso",
  "fofc",
  "two_fofc",
  "fmodel",
]

class residue_real_space(residue):
  # CC is 'score' attribute
  __slots__ = residue.__slots__ + __real_space_attr__

  @property
  def cc(self):
    return self.score

  @staticmethod
  def header():
    return "%-20s  %6s  %4s  %6s  %6s  %5s" % ("atom", "b_iso", "occ",
      "2Fo-Fc", "Fmodel", "CC")

  def as_string(self):
    return "%-20s  %6.2f  %4.2f  %6.2f  %6.2f  %5.3f" % (self.id_str(),
      self.b_iso, self.occupancy, self.two_fofc, self.fmodel, self.score)

  def as_table_row_phenix(self):
    return [ self.id_str(), self.b_iso, self.occupancy, self.two_fofc,
             self.fmodel, self.score ]

class data_statistics(slots_getstate_setstate):
  __slots__ = [
    "d_max",
    "d_min",
    "info",
    "wavelength",
    "n_refl",
    "n_refl_refine",
    "n_free",
    "r_work",
    "r_free",
    "twin_law",
    "wilson_b",
    "completeness",
    "r_work_outer",
    "r_free_outer",
    "d_max_outer",
    "d_min_outer",
    "completeness_outer",
    "n_refl_outer",
    "n_refl_refine_outer",
    "n_free_outer",
    "anomalous_flag",
  ]
  def __init__(self, fmodel, raw_data=None, n_bins=10,
      count_anomalous_pairs_separately=False):
    # FIXME n_bins should be automatic by default
    f_obs = fmodel.f_obs().deep_copy()
    r_free_flags = fmodel.r_free_flags().deep_copy()
    self.anomalous_flag = f_obs.anomalous_flag()
    self.d_max = f_obs.d_max_min()[0]
    self.d_min = f_obs.d_min()
    self.info = fmodel.info(n_bins=n_bins)
    self.wavelength = None
    if (raw_data is not None):
      self.wavelength = getattr(raw_data.info(), "wavelength", None)
      raw_data = raw_data.map_to_asu().eliminate_sys_absent()
      raw_data = raw_data.resolution_filter(d_max=self.d_max, d_min=self.d_min)
      if (raw_data.anomalous_flag() and not self.anomalous_flag):
        raw_data = raw_data.average_bijvoet_mates()
    else :
      raw_data = f_obs
    if (not count_anomalous_pairs_separately) and (self.anomalous_flag):
      f_obs = f_obs.average_bijvoet_mates()
      raw_data = raw_data.average_bijvoet_mates()
      r_free_flags = r_free_flags.average_bijvoet_mates()
    f_obs.setup_binner(n_bins=n_bins)
    self.n_refl = raw_data.indices().size()
    self.n_refl_refine = f_obs.indices().size()
    self.n_free = r_free_flags.data().count(True)
    self.r_free = self.info.r_free
    self.r_work = self.info.r_work
    self.twin_law = fmodel.twin_law
    self.wilson_b = fmodel.wilson_b()
    self.completeness = f_obs.completeness()
    # outer shell
    d_max_min_outer = f_obs.binner().bin_d_range(n_bins)
    assert (not None in d_max_min_outer)
    self.d_max_outer = d_max_min_outer[0]
    self.d_min_outer = d_max_min_outer[1]
    self.r_free_outer = fmodel.r_free(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    self.r_work_outer = fmodel.r_work(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    i_bin = f_obs.binner().range_used()[-1]
    outer_shell = f_obs.select(f_obs.binner().selection(i_bin))
    self.completeness_outer = outer_shell.completeness(d_max=self.d_max_outer)
    self.n_refl_outer = raw_data.resolution_filter(d_max=self.d_max_outer,
      d_min=self.d_min_outer).indices().size()
    f_obs_outer = f_obs.resolution_filter(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    self.n_refl_refine_outer = f_obs_outer.indices().size()
    r_free_flags_outer = r_free_flags.resolution_filter(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    self.n_free_outer = r_free_flags_outer.data().count(True)

  def show_summary(self, out=sys.stdout, prefix=""):
    print("%sHigh resolution       = %7.3f" % (prefix, self.d_min), file=out)
    print("%sR-work                = %8.4f" % (prefix, self.r_work), file=out)
    print("%sR-free                = %8.4f" % (prefix, self.r_free), file=out)

  def show(self, out=sys.stdout, prefix=""):
    def fv(fs, val):
      return format_value(fs, val, replace_none_with="----")
    if (not self.wavelength in [None, 0]):
      print("%sWavelength                 = %.4g" % (prefix,
        self.wavelength), file=out)
    print("%sResolution range           = %7.3f - %.3f (%.3f - %.3f)" \
      % (prefix, self.d_max, self.d_min, self.d_max_outer, self.d_min_outer), file=out)
    print("%sNumber of reflections      = %8d (%d)" % (prefix,
      self.n_refl, self.n_refl_outer), file=out)
    print("%s   after outlier rejection = %8d (%d)" % (prefix,
      self.n_refl_refine, self.n_refl_refine_outer), file=out)
    print("%sCompleteness               = %6.2f%% (%.2f%%)" % (prefix,
      self.info.completeness_in_range*100, self.completeness_outer*100), file=out)
    print("%sR-work                     = %8.4f (%s)" % (prefix,
      self.r_work, fv("%.4f", self.r_work_outer)), file=out)
    print("%sR-free                     = %8.4f (%s)" % (prefix,
      self.r_free, fv("%.4f", self.r_free_outer)), file=out)
    self.info.show_rwork_rfree_number_completeness(prefix=prefix, out=out,
      title="By resolution bin:")

class real_space(validation):
  """
  Real-space correlation calculation for residues
  """

  __slots__ = validation.__slots__ + \
              ['overall_rsc', 'fsc', 'everything', 'protein', 'other', 'water']
  program_description = "Analyze real space correlation"
  output_header = None
  gui_list_headers = [ "Residue", "B_iso", "Occupancy", "2Fo-Fc", "Fmodel", "CC" ]
  gui_formats = [ "%s", "%6.2f", "%4.2f", "%6.2f", "%6.2f", "%5.3f" ]
  wx_column_widths = [120] * 6

  def get_result_class(self) : return residue_real_space

  def __init__(self, model, fmodel, cc_min=0.8, molprobity_map_params=None):

    from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
    from mmtbx import real_space_correlation

    validation.__init__(self)

    pdb_hierarchy = model.get_hierarchy()
    crystal_symmetry = model.crystal_symmetry()

    # arrays for different components
    self.everything = list()
    self.protein = list()
    self.other = list()
    self.water = list()
    aa_codes = one_letter_given_three_letter

    # redo real_space_corelation.simple to use map objects instead of filenames
    self.overall_rsc = None
    rsc = None
    try :
      rsc_params = real_space_correlation.master_params().extract()
      rsc_params.detail="residue"
      rsc_params.map_1.fill_missing_reflections = False
      rsc_params.map_2.fill_missing_reflections = False
      use_maps = False
      if (molprobity_map_params is not None):
        rsc_params.map_coefficients_file_name = \
          molprobity_map_params.map_coefficients_file_name
        rsc_params.map_coefficients_label = \
          molprobity_map_params.map_coefficients_label
        if (molprobity_map_params.map_file_name is not None):
          use_maps = True
      # use mmtbx/command_line/map_model_cc.py for maps
      self.fsc = None
      if (use_maps):
        from iotbx.map_model_manager import map_model_manager
        from mmtbx.maps import map_model_cc
        #  XXX no longer exists from mmtbx.maps.import get_fsc
        from mmtbx.maps.mtriage import get_fsc  # XXX replaced above
        from iotbx.file_reader import any_file
        params = map_model_cc.master_params().extract()
        params.map_model_cc.resolution = molprobity_map_params.d_min
        map_object = any_file(molprobity_map_params.map_file_name).file_object

        # check that model crystal symmetry matches map crystal symmetry
        mmi = map_model_manager(
          map_manager    = map_object,
          model = model)

        rsc_object = map_model_cc.map_model_cc(
          mmi.map_data(), mmi.model().get_hierarchy(), mmi.crystal_symmetry(),
          params.map_model_cc)
        rsc_object.validate()
        rsc_object.run()
        rsc = rsc_object.get_results()
        self.overall_rsc = (rsc.cc_mask, rsc.cc_volume, rsc.cc_peaks)

        self.fsc = get_fsc(mmi.map_data(), mmi.model(), params.map_model_cc)
        self.fsc.atom_radius = rsc.atom_radius
        rsc = rsc.cc_per_residue
      # mmtbx/real_space_correlation.py for X-ray/neutron data and map
      # coefficients
      else:
        self.overall_rsc, rsc = real_space_correlation.simple(
          fmodel=fmodel,
          pdb_hierarchy=pdb_hierarchy,
          params=rsc_params,
          log=null_out())
    except Exception as e :
      raise
    else :
      assert ( (self.overall_rsc is not None) and (rsc is not None) )
      for i, result_ in enumerate(rsc):
        if (use_maps): # new rsc calculation (mmtbx/maps/model_map_cc.py)
          result = residue_real_space(
            chain_id=result_.chain_id,
            resname=result_.resname,
            resseq=result_.resseq,
            icode=result_.icode,
            altloc="",
            score=result_.cc,
            b_iso=result_.b_iso_mean,
            occupancy=result_.occ_mean,
            outlier=result_.cc < cc_min,
            xyz=result_.xyz_mean)
        else: # old rsc calculation (mmtbx/maps/real_space_correlation.py)
          result = residue_real_space(
            chain_id=result_.chain_id,
            resname=result_.residue.resname,
            resseq=result_.residue.resseq,
            icode=result_.residue.icode,
            altloc="",
            score=result_.cc,
            b_iso=result_.b,
            occupancy=result_.occupancy,
            fmodel=result_.map_value_1,
            two_fofc=result_.map_value_2,
            outlier=result_.cc < cc_min,
            xyz=result_.residue.atoms().extract_xyz().mean())
        if result.is_outlier():
          self.n_outliers += 1
        # XXX unlike other validation metrics, we always save the results for
        # the real-space correlation, since these are used as the basis for
        # the multi-criterion plot in Phenix.  The show() method will only
        # print outliers, however.
        if (result.resname != 'HOH'): # water is handled by waters.py
          self.everything.append(result)
          if result.resname in one_letter_given_three_letter:
            self.protein.append(result)
          else:
            self.other.append(result)
        self.everything += self.water
        self.results = self.protein

  def add_water(self, water=None):
    """
    Function for incorporating water results from water.py
    """
    if (water is not None):
      self.water = water
      self.everything += water

  def show_summary(self, out=sys.stdout, prefix=""):
    print(prefix +\
      "%d residues (including water) with CC(Fc,2mFo-DFc) < 0.8" % \
      self.n_outliers, file=out)

  def show(self, out=sys.stdout, prefix="  ", verbose=True):
    if (self.n_outliers > 0):
      print(prefix + self.get_result_class().header(), file=out)
      for result in (self.protein + self.other):
        if result.is_outlier():
          print(prefix + str(result), file=out)
    self.show_summary(out=out, prefix=prefix)

def merging_and_model_statistics(
    f_obs,
    f_model,
    r_free_flags,
    unmerged_i_obs,
    n_bins=20,
    sigma_filtering=Auto,
    anomalous=False,
    use_internal_variance=True):
  """
  Compute merging statistics - CC* in particular - together with measures of
  model quality using reciprocal-space data (R-factors and CCs).  See Karplus
  & Diederichs 2012 for rationale.
  """
  from iotbx import merging_statistics
  free_sel = r_free_flags
  # very important: must use original intensities for i_obs, not squared f_obs,
  # because French-Wilson treatment is one-way
  assert (unmerged_i_obs.sigmas() is not None)
  info = unmerged_i_obs.info()
  assert (info is not None)
  unmerged_i_obs = unmerged_i_obs.customized_copy(crystal_symmetry=f_obs)
  unmerged_i_obs = unmerged_i_obs.select(
    unmerged_i_obs.sigmas() >= 0).set_info(info)
  filter = merging_statistics.filter_intensities_by_sigma(
    array=unmerged_i_obs,
    sigma_filtering=sigma_filtering)
  i_obs = filter.array_merged
  unmerged_i_obs = filter.array
  # average Bijvoet pairs if not anomalous
  if (not anomalous):
    if (i_obs.anomalous_flag()):
      i_obs = i_obs.average_bijvoet_mates()
    if (f_obs.anomalous_flag()):
      f_obs = f_obs.average_bijvoet_mates()
    if (f_model.anomalous_flag()):
      f_model = f_model.average_bijvoet_mates()
    if (free_sel.anomalous_flag()):
      free_sel = free_sel.average_bijvoet_mates()
  # create Bijvoet pairs if an array is not anomalous
  else:
    if (not i_obs.anomalous_flag()):
      i_obs = i_obs.generate_bijvoet_mates()
    if (not f_obs.anomalous_flag()):
      f_obs = f_obs.generate_bijvoet_mates()
    if (not f_model.anomalous_flag()):
      f_model = f_model.generate_bijvoet_mates()
    if (not free_sel.anomalous_flag()):
      free_sel = free_sel.generate_bijvoet_mates()

  if (free_sel.data().count(True) == 0):
    raise Sorry("R-free array does not select any reflections.  To calculate "+
      "CC* and related statistics, a valid set of R-free flags must be used.")
  work_sel = free_sel.customized_copy(data=~free_sel.data())
  i_obs, f_obs = i_obs.common_sets(other=f_obs)
  i_obs, f_model = i_obs.common_sets(other=f_model)
  i_obs, work_sel = i_obs.common_sets(other=work_sel)
  i_obs, free_sel = i_obs.common_sets(other=free_sel)
  i_calc = abs(f_model).f_as_f_sq()
  d_max, d_min = i_calc.d_max_min()
  model_arrays = merging_statistics.model_based_arrays(
    f_obs=f_obs,
    i_obs=i_obs,
    i_calc=i_calc,
    work_sel=work_sel,
    free_sel=free_sel)
  return merging_statistics.dataset_statistics(
    i_obs=unmerged_i_obs,
    crystal_symmetry=i_calc,
    d_min=d_min,
    d_max=d_max,
    n_bins=n_bins,
    model_arrays=model_arrays,
    anomalous=anomalous,
    use_internal_variance=use_internal_variance,
    sigma_filtering=None) # no need, since it was done here
