
# XXX most of the functions in this module are deprectated and should be
# removed as soon as someone has time.

from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args
import libtbx.phil
import os
import sys
from six.moves import zip

#-----------------------------------------------------------------------
# MAP COEFFICIENT MANIPULATION

def create_map_from_pdb_and_mtz(
    pdb_file,
    mtz_file,
    output_file,
    fill=False,
    out=None,
    llg_map=False,
    remove_unknown_scatterering_type=False,
    assume_pdb_data=False):
  """
  Convenience function, used by phenix.fetch_pdb.

  :param remove_unknown_scatterering_type: if True, atoms with unknown
    scattering type will be removed before extracting the X-ray structure
    (which would crash otherwise).
  """
  if (out is None) : out = sys.stdout
  import iotbx.pdb
  from scitbx.array_family import flex
  pdb_in = iotbx.pdb.input(pdb_file)
  xrs = pdb_in.xray_structure_simple(enable_scattering_type_unknown=True)
  selection = xrs.scatterers().extract_scattering_types()!="unknown"
  if(selection.size()!=selection.count(True)):
    print("WARNING: removing %d atoms with unknown scattering type." %\
      selection.count(False), file=out)
  xrs = xrs.select(selection)
  fast_maps_from_hkl_file(
    file_name=mtz_file,
    xray_structure=xrs,
    map_out=output_file,
    log=out,
    auto_run=True,
    quiet=True,
    anomalous_map=True,
    llg_map=llg_map,
    fill_maps=fill,
    assume_pdb_data=assume_pdb_data)

class fast_maps_from_hkl_file(object):
  """
  Automation wrapper for calculating standard 2mFo-DFc, mFo-DFc, and optional
  anomalous difference map coefficients, starting from an X-ray structure and
  an MTZ file name.  This will attempt to guess the correct data array
  automatically (giving preference to anomalous amplitudes).

  The organization of this class is somewhat messy, as it is designed to be
  run in parallel for multiple structures.
  """
  def __init__(self,
                file_name,
                xray_structure,
                scattering_table="wk1995",
                f_label=None,
                r_free_label=None,
                map_out=None,
                log=sys.stdout,
                auto_run=True,
                quiet=False,
                anomalous_map=False,
                fill_maps=True,
                save_fmodel=False,
                llg_map=False,
                assume_pdb_data=False, # FIXME should be default, needs testing
                ):
    adopt_init_args(self, locals())
    from iotbx import file_reader
    from scitbx.array_family import flex
    if f_label is None and not quiet:
      print("      no label for %s, will try default labels" % \
        os.path.basename(file_name), file=log)
    f_obs = r_free_flags = None
    # FIXME this whole block needs to disappear
    if (not assume_pdb_data):
      fallback_f_obs = []
      default_labels = ["F(+),SIGF(+),F(-),SIGF(-)",
        "I(+),SIGI(+),I(-),SIGI(-)",
        "F,SIGF","FOBS,SIGFOBS", "FOBS_X", "IOBS,SIGIOBS"]
      default_rfree_labels = ["FreeR_flag", "FREE", "R-free-flags"]
      all_labels = []
      best_label = sys.maxsize
      data_file = file_reader.any_file(file_name, force_type="hkl")
      for miller_array in data_file.file_server.miller_arrays :
        labels = miller_array.info().label_string()
        all_labels.append(labels)
        if (labels == f_label):
          f_obs = miller_array
          break
        elif (f_label is None) and (labels in default_labels):
          label_score = default_labels.index(labels)
          if (label_score < best_label):
            f_obs = miller_array
            best_label = label_score
        elif miller_array.is_xray_amplitude_array():
          fallback_f_obs.append(miller_array)
      if f_obs is None :
        if (len(fallback_f_obs) == 1) and (f_label is None):
          for array in fallback_f_obs :
            if (array.anomalous_flag()) and (self.anomalous_map):
              f_obs = array
              break
          else :
            f_obs = fallback_f_obs[0]
        else :
          raise Sorry(("Couldn't find %s in %s.  Please specify valid "+
            "column labels (possible choices: %s)") % (f_label, file_name,
              " ".join(all_labels)))
      if (f_obs.is_xray_intensity_array()):
        f_obs = f_obs.f_sq_as_f()
      sys_abs_flags = f_obs.sys_absent_flags().data()
      f_obs = f_obs.map_to_asu().select(selection=~sys_abs_flags)
      r_free = data_file.file_server.get_r_free_flags(
        file_name=None,
        label=r_free_label,
        test_flag_value=None,
        parameter_scope=None,
        disable_suitability_test=False,
        return_all_valid_arrays=True)
      if (len(r_free) == 0):
        self.f_obs = f_obs
        self.r_free_flags = f_obs.array(
          data=flex.bool(f_obs.data().size(),False))
      else :
        array, test_flag_value = r_free[0]
        new_flags = array.customized_copy(
          data=array.data() == test_flag_value).map_to_asu()
        if (f_obs.anomalous_flag()) and (not new_flags.anomalous_flag()):
          new_flags = new_flags.generate_bijvoet_mates()
        self.r_free_flags = new_flags.common_set(f_obs)
        self.f_obs = f_obs.common_set(self.r_free_flags)
    else :
      import mmtbx.wwpdb.utils
      data_info = mmtbx.wwpdb.utils.find_data_arrays(
        mtz_file=file_name, log=null_out())
      self.f_obs, self.r_free_flags = data_info.data_and_flags()
    self.log = None
    self.fmodel = None
    if auto_run :
      self.run()

  def get_maps_from_fmodel(self):
    import mmtbx.utils
    from scitbx.array_family import flex
    mmtbx.utils.setup_scattering_dictionaries(
      scattering_table = "n_gaussian",
      xray_structure   = self.xray_structure,
      d_min            = self.f_obs.d_min(),
      log              = null_out())
    fmodel = mmtbx.utils.fmodel_simple(
      xray_structures=[self.xray_structure],
      scattering_table = self.scattering_table,
      f_obs=self.f_obs,
      r_free_flags=self.r_free_flags,
      outliers_rejection=True,
      skip_twin_detection=False,
      bulk_solvent_correction=True,
      anisotropic_scaling=True)
    if (self.save_fmodel):
      self.fmodel = fmodel
    (f_map, df_map) = get_maps_from_fmodel(fmodel)
    anom_map = None
    if (self.anomalous_map) and (self.f_obs.anomalous_flag()):
      anom_map = get_anomalous_map(fmodel)
    return f_map, df_map, anom_map

  def run(self):
    import iotbx.map_tools
    (f_map, df_map, anom_map) = self.get_maps_from_fmodel()
    if self.map_out is None :
      self.map_out = os.path.splitext(self.file_name)[0] + "_map_coeffs.mtz"
    iotbx.map_tools.write_map_coeffs(f_map, df_map, self.map_out, anom_map)

def get_maps_from_fmodel(fmodel, fill_missing_f_obs=True,
    exclude_free_r_reflections=False):
  map_manager = fmodel.electron_density_map()
  two_fofc_coeffs = map_manager.map_coefficients(map_type = "2mFo-DFc",
    exclude_free_r_reflections=exclude_free_r_reflections,
    fill_missing=fill_missing_f_obs)
  if two_fofc_coeffs.anomalous_flag():
    two_fofc_coeffs = two_fofc_coeffs.average_bijvoet_mates()
  fofc_coeffs = map_manager.map_coefficients(map_type = "mFo-DFc",
    exclude_free_r_reflections=exclude_free_r_reflections,
    fill_missing=fill_missing_f_obs)
  if fofc_coeffs.anomalous_flag():
    fofc_coeffs = fofc_coeffs.average_bijvoet_mates()
  return (two_fofc_coeffs, fofc_coeffs)

def get_anomalous_map(fmodel):
  map_manager = fmodel.electron_density_map()
  anom_coeffs = map_manager.map_coefficients(map_type="anom")
  if (anom_coeffs.anomalous_flag()):
    anom_coeffs = anom_coeffs.average_bijvoet_mates()
  return anom_coeffs

class generate_water_omit_map(object):
  """
  Calculate standard map coefficients (with R-free flags omitted, and missing
  reflections filled) with water occupancies set to zero.  Used for ligand
  fitting after refinement.
  """
  def __init__(self,
      fmodel,
      pdb_hierarchy,
      skip_if_no_waters=False,
      exclude_free_r_reflections=False,
      fill_missing_f_obs=False,
      write_f_model=False,
      log=None):
    if (log is None):
      log = null_out()
    sel_cache = pdb_hierarchy.atom_selection_cache()
    water_sel = sel_cache.selection("resname HOH")
    xrs = fmodel.xray_structure
    self.crystal_symmetry = xrs.crystal_symmetry()
    assert (water_sel.size() == xrs.scatterers().size())
    self.n_waters = water_sel.count(True)
    self.two_fofc_map = self.fofc_map = self.anom_map = self.fmodel_map = None
    if (self.n_waters == 0) and (skip_if_no_waters):
      print("  No waters in model, skipping omit map calculation", file=log)
      self.pdb_hierarchy_omit = pdb_hierarchy
    else :
      print("  Calculating new 2mFo-DFc and mFo-DFc maps...", file=log)
      print("    %d waters temporarily set to zero occupancy" % \
        self.n_waters, file=log)
      print("", file=log)
      self.pdb_hierarchy_omit = pdb_hierarchy.select(~water_sel)
      occ = xrs.scatterers().extract_occupancies()
      if (self.n_waters > 0):
        xrs.set_occupancies(value=0, selection=water_sel)
      fmodel.update_xray_structure(xrs,
        update_f_calc=True)
    #    update_f_mask=False)
      fmodel.update_all_scales(log=log)
      self.two_fofc_map = fmodel.electron_density_map(
        ).map_coefficients(
          map_type="2mFo-DFc",
          exclude_free_r_reflections=exclude_free_r_reflections,
          fill_missing=fill_missing_f_obs)
      self.fofc_map = fmodel.electron_density_map(
        ).map_coefficients(
          map_type="mFo-DFc",
          exclude_free_r_reflections=exclude_free_r_reflections,
          fill_missing=False,
          merge_anomalous=True)
      if (fmodel.f_obs().anomalous_flag()):
        self.anom_map = fmodel.electron_density_map(
          ).map_coefficients(
            map_type="anom",
            exclude_free_r_reflections=exclude_free_r_reflections,
            fill_missing=fill_missing_f_obs)
      if (write_f_model):
        self.fmodel_map = fmodel.f_model().average_bijvoet_mates()

  def write_map_coeffs(self, file_name):
    assert (self.two_fofc_map is not None)
    write_map_coeffs(
      fwt_coeffs=self.two_fofc_map,
      delfwt_coeffs=self.fofc_map,
      anom_coeffs=self.anom_map,
      fmodel_coeffs=self.fmodel_map,
      file_name=file_name)

  def write_pdb_file(self, file_name):
    self.pdb_hierarchy_omit.write_pdb_file(
      file_name=file_name,
      crystal_symmetry=self.crystal_symmetry)

def get_default_map_coefficients_phil(
    show_filled=True,
    disable_filled=False,
    exclude_free_r_reflections=False,
    enclosing_scope=None,
    as_phil_object=True):
  """
  Generate the map coefficients PHIL string for phenix.maps, phenix.refine,
  and related applications, based on several toggles which control the
  handling of missing reflections and R-free flags.  Used by the Phenix GUI.
  """
  secondary_two_fofc_params = ""
  if (not disable_filled):
    if (show_filled):
      secondary_two_fofc_params = """
      map_coefficients {
        mtz_label_amplitudes = 2FOFCWT_no_fill
        mtz_label_phases = PH2FOFCWT_no_fill
        map_type = 2mFo-DFc
        fill_missing_f_obs = False
        exclude_free_r_reflections = %s
      }""" % exclude_free_r_reflections
    else :
      secondary_two_fofc_params = """
      map_coefficients {
        mtz_label_amplitudes = 2FOFCWT_filled
        mtz_label_phases = PH2FOFCWT_filled
        map_type = 2mFo-DFc
        fill_missing_f_obs = True
        exclude_free_r_reflections = %s
      }""" % exclude_free_r_reflections
  maps_phil = """
    map_coefficients {
      mtz_label_amplitudes = 2FOFCWT
      mtz_label_phases = PH2FOFCWT
      map_type = 2mFo-DFc
      fill_missing_f_obs = %(show_filled)s
      exclude_free_r_reflections = %(exclude_rfree)s
    }
    %(other_two_fofc)s
    map_coefficients {
      mtz_label_amplitudes = FOFCWT
      mtz_label_phases = PHFOFCWT
      map_type = mFo-DFc
      exclude_free_r_reflections = %(exclude_rfree)s
    }
    map_coefficients {
      map_type = anomalous
      format = *mtz phs
      mtz_label_amplitudes = ANOM
      mtz_label_phases = PHANOM
      exclude_free_r_reflections = %(exclude_rfree)s
    }
""" % { "exclude_rfree" : exclude_free_r_reflections,
        "other_two_fofc" : secondary_two_fofc_params,
        "show_filled" : show_filled and not disable_filled, }
  if (enclosing_scope is not None):
    maps_phil = """%s {\n%s\n}""" % (enclosing_scope, maps_phil)
  if (as_phil_object):
    return libtbx.phil.parse(maps_phil)
  return maps_phil

# XXX redundant, needs to be eliminated
def write_map_coeffs(*args, **kwds):
  import iotbx.map_tools
  return iotbx.map_tools.write_map_coeffs(*args, **kwds)

#-----------------------------------------------------------------------
# XPLOR MAP OUTPUT

# TODO: make more modular!
def write_xplor_map_file(coeffs, frac_min, frac_max, file_base):
  fft_map = coeffs.fft_map(resolution_factor=1/3.0)
  fft_map.apply_sigma_scaling()
  n_real = fft_map.n_real()
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  title_lines=["REMARK map covering model + 3.0A buffer"]
  file_name = "%s.map" % file_base
  fft_map.as_xplor_map(
    file_name=file_name,
    title_lines=title_lines,
    gridding_first=gridding_first,
    gridding_last=gridding_last)
  return file_name

def write_xplor_map(sites_cart, unit_cell, map_data, n_real, file_name,
    buffer=10):
  import iotbx.xplor.map
  if sites_cart is not None :
    frac_min, frac_max = unit_cell.box_frac_around_sites(
      sites_cart=sites_cart,
      buffer=buffer)
  else :
    frac_min, frac_max = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  gridding = iotbx.xplor.map.gridding(n     = map_data.focus(),
                                      first = gridding_first,
                                      last  = gridding_last)
  iotbx.xplor.map.writer(
    file_name          = file_name,
    is_p1_cell         = True,
    title_lines        = [' None',],
    unit_cell          = unit_cell,
    gridding           = gridding,
    data               = map_data,
    average            = -1,
    standard_deviation = -1)

# XXX backwards compatibility
# TODO remove these ASAP once GUI is thoroughly tested
def extract_map_coeffs(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.extract_map_coeffs(*args, **kwds)

def map_coeffs_from_mtz_file(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.map_coeffs_from_mtz_file(*args, **kwds)

def extract_phenix_refine_map_coeffs(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.extract_phenix_refine_map_coeffs(*args,
    **kwds)

def get_map_coeff_labels(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.get_map_coeff_labels(*args, **kwds)

def get_map_coeffs_for_build(server):
  return get_map_coeff_labels(server, build_only=True)

def format_map_coeffs_for_resolve(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.format_map_coeffs_for_resolve(*args,
    **kwds)

def decode_resolve_map_coeffs(*args, **kwds):
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.decode_resolve_map_coeffs(*args, **kwds)
