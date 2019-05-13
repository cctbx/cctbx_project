
from __future__ import division
from __future__ import print_function
from libtbx import Auto, adopt_init_args
from libtbx import easy_mp
import os.path
import sys

master_phil_str = """
selection = all
  .type = atom_selection
map_type = mFo-DFc
  .type = str
occ = 0.5
  .type = float
remove_waters = True
  .type = bool
omit_fraction = 0.02
  .type = float
optimize_binning = True
  .type = bool
box_cushion_radius = 2.5
  .type = float
fill_missing_f_obs = False
  .type = bool
exclude_free_r_reflections = True
  .type = bool
resolution_factor = 0.25
  .type = float
flatten_background = False
  .type = bool
nproc = Auto
  .type = int
output {
  mtz_file = partial_omit_coeffs.mtz
    .type = path
  ccp4_map = None
    .type = path
  verbose = False
    .type = bool
}
"""

map_type_labels = {
  "2mFo-DFc" : "2FOFCWT",
  "mFo-DFc" : "FOFCWT",
  "anom" : "ANOM",
  "anom_residual" : "ANOM_DIFF",
  "llg" : "LLG",
}

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string=master_phil_str,
    enable_automatic_twin_detection=True)

class partial_omit_map(object):
  def __init__(self,
      fmodel,
      selection,
      occupancy,
      map_type="mFo-DFc",
      omit_fraction=0.02,
      selection_delete=None,
      fill_missing_f_obs=False,
      exclude_free_r_reflections=False,
      optimize_binning=True,
      box_cushion_radius=2.5,
      nproc=Auto,
      out=sys.stdout):
    adopt_init_args(self, locals())
    import mmtbx.maps.composite_omit_map
    occ_saved = fmodel.xray_structure.scatterers().extract_occupancies()
    if (omit_fraction == 1.0):
      self.omit_groups = [ mmtbx.maps.composite_omit_map.omit_regions(
        serial=1,
        selection=selection) ]
    else :
      self.omit_groups = mmtbx.maps.composite_omit_map.create_omit_regions(
        xray_structure=fmodel.xray_structure,
        selection=selection,
        fraction_omit=self.omit_fraction,
        optimize_binning=optimize_binning,
        box_cushion_radius=box_cushion_radius,
        log=out)
      for group in self.omit_groups :
        group.show(out=out)
    self.omit_map_coeffs = easy_mp.pool_map(
      fixed_func=self,
      iterable=self.omit_groups,
      processes=nproc)
    fmodel.xray_structure.scatterers().set_occupancies(occ_saved)
    fmodel.update_xray_structure(update_f_calc=True)

  def combine_maps(self,
      flatten_background=False,
      sigma_scaling=False,
      resolution_factor=0.25):
    if (self.omit_fraction == 1.0):
      self.map_coeffs = self.omit_map_coeffs[0]
      fft_map = self.map_coeffs.fft_map(
        resolution_factor=resolution_factor).apply_volume_scaling()
      if (sigma_scaling):
        fft_map.apply_sigma_scaling()
      self.composite_map = fft_map.real_map_unpadded()
    else :
      import mmtbx.maps.composite_omit_map
      background_map_coeffs = self.fmodel.map_coefficients(
        map_type=self.map_type,
        fill_missing=self.fill_missing_f_obs,
        exclude_free_r_reflections=self.exclude_free_r_reflections,
        merge_anomalous=True)
      self.composite_map = mmtbx.maps.composite_omit_map.combine_maps(
        map_arrays=self.omit_map_coeffs,
        omit_groups=self.omit_groups,
        background_map_coeffs=background_map_coeffs,
        flatten_background=flatten_background,
        resolution_factor=resolution_factor,
        sigma_scaling=sigma_scaling)
      composite_map_coeffs = background_map_coeffs.structure_factors_from_map(
        map=self.composite_map,
        use_sg=True)
      self.map_coeffs = composite_map_coeffs
    return self.map_coeffs

  def __call__(self, omit_group):
    from scitbx.array_family import flex
    fmodel_tmp = self.fmodel.deep_copy()
    xrs_omit = fmodel_tmp.xray_structure.deep_copy_scatterers()
    xrs_omit.set_occupancies(self.occupancy, selection=omit_group.selection)
    if (self.selection_delete is not None):
      selection_partial = flex.bool(xrs_omit.scatterers().size(), False)
      selection_partial.set_selected(omit_group.selection, True)
      selection_delete = selection_partial & self.selection_delete
      n_del = selection_delete.count(True)
      if (n_del > 0):
        #print "setting %d waters to zero occupancy" % n_del
        xrs_omit.set_occupancies(0., selection=selection_delete)
    fmodel_tmp.update_xray_structure(xrs_omit,
      update_f_mask=False,
      update_f_calc=True)
    return fmodel_tmp.map_coefficients(
      map_type=self.map_type,
      fill_missing=self.fill_missing_f_obs,
      exclude_free_r_reflections=self.exclude_free_r_reflections,
      merge_anomalous=True)

def run(args, out=sys.stdout):
  import mmtbx.command_line
  import iotbx.map_tools
  from scitbx.array_family import flex
  usage_string = """\
mmtbx.partial_omit_map model.pdb data.mtz [occ=0.5] [omit_fraction=0.02] [...]

Generate a simple (un-refined) composite omit map with atoms at reduced rather
than zero occupancy.  Intended for use in the analysis of static disorder -
note that no treatment of phase bias is performed.
"""
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=False,
    prefer_anomalous=True,
    usage_string=usage_string,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    out=out)
  fmodel = cmdline.fmodel
  xray_structure = fmodel.xray_structure
  pdb_hierarchy = cmdline.pdb_hierarchy
  params = cmdline.params
  sel_cache = pdb_hierarchy.atom_selection_cache()
  selection = sel_cache.selection(params.selection)
  assert (selection.count(True) > 0)
  selection_delete = flex.bool(selection.size(), False)
  if (params.remove_waters):
    selection_delete = sel_cache.selection("resname HOH")
  map_driver = partial_omit_map(
    fmodel=fmodel,
    selection=selection,
    selection_delete=selection_delete,
    occupancy=params.occ,
    map_type=params.map_type,
    omit_fraction=params.omit_fraction,
    fill_missing_f_obs=params.fill_missing_f_obs,
    exclude_free_r_reflections=params.exclude_free_r_reflections,
    optimize_binning=params.optimize_binning,
    box_cushion_radius=params.box_cushion_radius,
    nproc=params.nproc,
    out=out)
  if (params.output.mtz_file is not None):
    map_coeffs = map_driver.combine_maps(
      flatten_background=params.flatten_background,
      sigma_scaling=False,
      resolution_factor=params.resolution_factor)
    iotbx.map_tools.write_map_coefficients_generic(
      map_coeffs=[map_coeffs],
      map_types=[params.map_type],
      file_name=params.output.mtz_file)
    print("Wrote map coefficients to %s" % params.output.mtz_file, file=out)
  if (params.output.ccp4_map is not None):
    map_driver.combine_maps(
      flatten_background=params.flatten_background,
      sigma_scaling=True,
      resolution_factor=params.resolution_factor)
    iotbx.map_tools.write_ccp4_map(
      sites_cart=fmodel.xray_structure.sites_cart(),
      unit_cell=fmodel.xray_structure.unit_cell(),
      map_data=map_driver.composite_map,
      n_real=map_driver.composite_map.focus(),
      file_name=params.output.ccp4_map)
    print("Wrote CCP4 map to %s" % params.output.ccp4_map, file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
