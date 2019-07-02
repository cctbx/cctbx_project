
"""
For the obsessive-compulsive who absolutely must have perfectly contoured
density around their ligand.  Handy for making figures, but prone to abuse,
which is why it's not in the command_line directory.
"""

from __future__ import absolute_import, division, print_function
import sys

def flatten_map(map, xray_structure, selection):
  from cctbx import maptbx
  from scitbx.array_family import flex
  sites = xray_structure.sites_cart().select(selection)
  hd_sel = xray_structure.hd_selection()
  radii = flex.double()
  for i_seq in selection :
    if (hd_sel[i_seq]):
      radii.append(1.0)
    else :
      radii.append(1.5)
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = xray_structure.unit_cell(),
    fft_n_real = map.focus(),
    fft_m_real = map.all(),
    sites_cart = sites,
    site_radii = radii)
  bg_sel = flex.bool(map.size(), True)
  bg_sel.set_selected(sel, False)
  map.as_1d().set_selected(bg_sel, 0)
  return map

def run(args, out=sys.stdout):
  master_phil_str = """
map_coeffs = None
  .type = path
model = None
  .type = path
selection = all
  .type = atom_selection
prefix = flattened
  .type = str
write_pymol = False
  .type = bool
"""
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    reflection_file_def="map_coeffs",
    pdb_file_def="model",
    usage_string="""\
phenix.flatten_map_outside_selection model.pdb maps.mtz selection='chain A'

For each set of map coefficients in the input MTZ file, performs an FFT to
obtain the electron density, and sets all grid points outside the defined
atom selection to zero.  Used for generating figures in PyMOL.""")
  params = cmdline.work.extract()
  assert (not None in [params.model, params.map_coeffs, params.selection])
  from iotbx import file_reader
  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  sites_cart = xrs.sites_cart()
  sel = hierarchy.atom_selection_cache().selection(params.selection)
  isel = sel.iselection()
  mtz_in = file_reader.any_file(params.map_coeffs, force_type="hkl")
  two_fofc = fofc = None
  for miller_array in mtz_in.file_server.miller_arrays :
    if (miller_array.is_complex_array()):
      label = miller_array.info().labels[0]
      if (label.startswith("F-model")) : continue
      real_map = miller_array.fft_map(resolution_factor=0.25
        ).apply_sigma_scaling().real_map_unpadded()
      flatten_map(
        map=real_map,
        xray_structure=xrs,
        selection=isel)
      import iotbx.map_tools
      file_name = params.prefix + "_%s.ccp4" % label
      iotbx.map_tools.write_ccp4_map(
        sites_cart=sites_cart,
        unit_cell=xrs.unit_cell(),
        map_data=real_map,
        n_real=real_map.focus(),
        file_name=file_name,
        buffer=5)
      print("wrote %s" % file_name)
      if (label == "2FOFCWT"):
        two_fofc = file_name
      elif (label == "FOFCWT"):
        fofc = file_name
  if (params.write_pymol):
    f = open(params.prefix + ".pml", "w")
    f.write("""set normalize_ccp4_maps, 0\n""")
    f.write("""load %s\n""" % params.model)
    f.write("""show sticks, (%s) and not elem H\n""" % params.selection)
    f.write("""color yellow, elem C\n""")
    f.write("""hide lines\n""")
    f.write("""hide nonbonded\n""")
    f.write("""zoom %s\n""" % params.selection)
    if (two_fofc is not None):
      f.write("""load %s, 2fofc\n""" % two_fofc)
      f.write("""isomesh m1, 2fofc, 1.0, %s, 5\n""" % params.selection)
      f.write("""color grey50, m1\n""")
    if (fofc is not None):
      f.write("""load %s, fofc\n""" % fofc)
      f.write("""isomesh m2, fofc, 3.0, %s, 5\n""" % params.selection)
      f.write("""color green, m2\n""")
    f.close()

if (__name__ == "__main__"):
  run(sys.argv[1:])
