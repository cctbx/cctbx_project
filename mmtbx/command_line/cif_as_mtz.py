# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_mtz

import sys, os, re
import string
from cctbx.array_family import flex
from libtbx import runtime_utils
from cctbx import crystal
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import crystal_symmetry_from_pdb
import iotbx.mtz
from libtbx import smart_open
import mmtbx.utils
import iotbx.pdb
import iotbx.phil

"""
Notes on CIF (source: http://www.ccp4.ac.uk/html/mtz2various.html)
All reflections in the MTZ input file will be output to the CIF file. However,
there are ways to flag certain reflections with the data type _refln.status.
Observed reflections will be flagged with 'o'. Unobserved reflections, i.e.
those flagged as missing in the relevant amplitude or intensity column, will be
flagged as 'x'; these reflections will not be added to _reflns.number_obs. The
'free' reflections will be flagged as 'f'. The keyword FREEVAL can be used to
indicate this set. Systematically absent reflections are flagged with '-'.

If the RESO keyword is specified then reflections at higher or lower resolution
than the limits given, will be written with _refln.status 'h' or 'l'
respectively. The limits will be written to the CIF as the values of
_refine.ls_d_res_high and _refine.ls_d_res_low.

If EXCLUDE SIG is given then reflections for which F < <value>*sigma(F), and
which satisfy the resolution limits (if given), will be written with
_refln.status '<'. The value of _reflns.number_obs excludes all reflections
which do not satisfy the condition on sigma(F). All other sub-keywords of
EXCLUDE are ignored for CIF output.
NB: The translation of the RESOLUTION and EXCLUDE SIGP conditions to
_refln.status values does not imply that the the use of these conditions is
good crystallographic practice. Be prepared to justify why you have excluded
any data from your final refinement!
"""

def run(args, command_name = "phenix.cif_as_mtz"):
  if (len(args) == 0): args = ["--help"]
  try:
    command_line = (iotbx_option_parser(
      usage="%s [reflection_cif_file] [options]" % command_name,
      description='Example: %s r1o9ksf.ent --symmetry=pdb1o9k.ent'%command_name)
      .enable_symmetry_comprehensive()
      .option(None, "--output_file_name",
        action="store",
        default=False,
        type="string",
        help="Output mtz file name.")
      .option(None, "--use_model",
        action="store",
        default=False,
        type="string",
        help="Use PDB model to make better guess about reflection data type.")
      .option(None, "--wavelength_id",
        action="store",
        default=None,
        type="int",
        help="Extract data set with given wavelength_id.")
      .option(None, "--crystal_id",
        action="store",
        default=None,
        type="int",
        help="Extract data set with given crystal_id.")
      .option(None, "--output_r_free_label",
        action="store",
        default="R-free-flags",
        type="string",
        help="MTZ column label to use for R-free flags (default: R-free-flags)")
      .option(None, "--merge",
        action="store_true",
        help="Merge non-unique data where present.")
      .option(None, "--remove_systematic_absences",
        action="store_true",
        help="Remove systematic absent reflections.")
      .option(None, "--map_to_asu",
        action="store_true",
        help="Map to asymmetric unit.")
      .option("--show_details_if_error",
          action="store_true",
          help="Show data details for some errors.")
      .option("--show_log",
          action="store_true",
          help="Show some output.")
    ).process(args=args)
  except Exception, e:
    if(str(e) != "0"): print str(e)
    sys.exit(0)
  crystal_symmetry = command_line.symmetry
  if(command_line.symmetry.unit_cell() is None or
     command_line.symmetry.space_group_info() is None):
    if(command_line.options.use_model):
      crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
         file_name=command_line.options.use_model)
  if(crystal_symmetry.unit_cell() is None or
     crystal_symmetry.space_group_info() is None):
    raise Sorry(
      "Crystal symmetry is not defined. Please use the --symmetry option.\n"
      "Type %s without arguments to see more options."%command_name)
  if(len(command_line.args) > 1):
    print "%d arguments are given from the command line:"% \
      len(command_line.args), command_line.args
    raise Sorry("Please specify one reflection cif file.")
  file_name = command_line.args[0]
  if(not os.path.isfile(file_name)):
    raise Sorry("File is not found: %s"%file_name)
  output_r_free_label = command_line.options.output_r_free_label
  if ((not output_r_free_label[0] in string.uppercase) or
      (re.search("[^a-zA-Z0-9_\-]", output_r_free_label))) :
    raise Sorry(("%s is not a suitable column label.  MTZ format requires "+
      "an uppercase letter as the first character, and only alphanumeric "+
      "characters or hyphens in the rest of the string.")% output_r_free_label)
  process_files(
    file_name=file_name,
    crystal_symmetry=crystal_symmetry,
    pdb_file_name=command_line.options.use_model,
    output_file_name=command_line.options.output_file_name,
    wavelength_id=command_line.options.wavelength_id,
    crystal_id=command_line.options.crystal_id,
    show_details_if_error=command_line.options.show_details_if_error,
    output_r_free_label=command_line.options.output_r_free_label,
    merge_non_unique_under_symmetry=command_line.options.merge,
    map_to_asu=command_line.options.map_to_asu,
    remove_systematic_absences=command_line.options.remove_systematic_absences)

def process_files (file_name,
                   crystal_symmetry,
                   pdb_file_name,
                   output_file_name,
                   wavelength_id,
                   crystal_id,
                   show_details_if_error,
                   output_r_free_label,
                   merge_non_unique_under_symmetry=False,
                   map_to_asu=False,
                   remove_systematic_absences=False) :
  mtz_object = extract(
    file_name                       = file_name,
    crystal_symmetry                = crystal_symmetry,
    wavelength_id                   = wavelength_id,
    crystal_id                      = crystal_id,
    show_details_if_error           = show_details_if_error,
    output_r_free_label             = output_r_free_label,
    merge_non_unique_under_symmetry = merge_non_unique_under_symmetry,
    map_to_asu                      = map_to_asu,
    remove_systematic_absences      = remove_systematic_absences)
  if(mtz_object is not None):
    if (pdb_file_name):
      pdb_raw_records = smart_open.for_reading(
        file_name=pdb_file_name).read().splitlines()
      xray_structure = None
      try:
        xray_structure = iotbx.pdb.input(file_name =
          pdb_file_name).xray_structure_simple()
      except Exception, e:
        print "Cannot extract xray_structure: ", str(e)
      if(xray_structure is not None):
        miller_arrays = mtz_object.as_miller_arrays()
        if(len(miller_arrays) == 1):
          f_obs = miller_arrays[0]
          r_free_flags = None
        else:
          if (len(miller_arrays) != 2) :
            raise Sorry("The --use-model flag is not supported when more "+
              "than one array of experimental data is present.")
          r_free_flags = None
          for miller_array in miller_arrays:
            if(miller_array.observation_type() is not None):
              f_obs = miller_array
            else:
              r_free_flags = miller_array
              assert isinstance(r_free_flags.data(), flex.int)
        data_label = f_obs.info().labels[0]
        if(r_free_flags is not None):
          f_obs = f_obs.common_set(r_free_flags)
          r_free_flags = r_free_flags.common_set(f_obs)
        f_obs = f_obs.merge_equivalents().array()
        if(r_free_flags is not None):
          r_free_flags = r_free_flags.merge_equivalents().array()
        mtz_object = mmtbx.utils.guess_observation_type(
          f_obs          = f_obs,
          label          = data_label,
          xray_structure = xray_structure,
          r_free_flags   = r_free_flags).mtz_object()
    if not output_file_name :
      basename = os.path.basename(file_name)
      if(basename[-4:-3] == "."): output_file_name = basename[:-4]+".mtz"
      elif(basename[-5:-4] == "."): output_file_name = basename[:-5]+".mtz"
      elif(basename.endswith(".ent.gz")): output_file_name=basename[:-7]+".mtz"
      else: output_file_name = basename+".mtz"
    mtz_object.write(file_name = output_file_name)
    return mtz_object.n_reflections()

def extract(file_name,
            crystal_symmetry,
            wavelength_id,
            crystal_id,
            show_details_if_error,
            output_r_free_label,
            merge_non_unique_under_symmetry,
            map_to_asu,
            remove_systematic_absences):
  import iotbx.cif
  all_miller_arrays = iotbx.cif.reader(file_path=file_name).build_miller_arrays()
  if (len(all_miller_arrays) == 0) :
    raise Sorry("No data arrays were found in this CIF file.  Please make "+
      "sure that the file contains reflection data, rather than the refined "+
      "model.")
  column_labels = set()

  def get_label(miller_array):
    label = None
    for l in miller_array.info().labels:
      if ('_meas' in l) :
        if miller_array.is_xray_amplitude_array():
          label = "FOBS"
        elif miller_array.is_xray_intensity_array():
          label = "IOBS"
        break
      elif miller_array.anomalous_flag():
        if miller_array.is_xray_amplitude_array():
          label = "F"
        elif miller_array.is_xray_intensity_array():
          label = "I"
        break
      elif 'status' in l or '_free' in l:
        label = output_r_free_label
        break
      elif miller_array.is_hendrickson_lattman_array():
        label = "HL"
    if label is not None:
      label_base = label
      i = 1
      while label in column_labels:
        label = label_base + "-%i" %(i)
        i += 1
    return label

  mtz_object = iotbx.mtz.object() \
    .set_title(title="phenix.cif_as_mtz") \
    .set_space_group_info(space_group_info=crystal_symmetry.space_group_info())
  unit_cell=crystal_symmetry.unit_cell()
  mtz_crystals = {}
  mtz_object.set_hkl_base(unit_cell=unit_cell)
  from iotbx.reflection_file_utils import cif_status_flags_as_int_r_free_flags
  for i, (data_name, miller_arrays) in enumerate(all_miller_arrays.iteritems()):
    for ma in miller_arrays.values():
      other_symmetry = crystal_symmetry
      try:
        crystal_symmetry = other_symmetry.join_symmetry(
          other_symmetry=ma.crystal_symmetry(),
          force=True)
      except AssertionError, e:
        str_e = str(e)
        from cStringIO import StringIO
        s = StringIO()
        if "Space group is incompatible with unit cell parameters." in str_e:
          other_symmetry.show_summary(f=s)
          ma.crystal_symmetry().show_summary(f=s)
          str_e += "\n%s" %(s.getvalue())
          raise Sorry(str_e)
        else:
          raise
      ma = ma.customized_copy(
        crystal_symmetry=crystal_symmetry).set_info(ma.info())
      labels = ma.info().labels
      label = get_label(ma)
      if label is None: continue
      elif label.startswith(output_r_free_label):
        ma, _ = cif_status_flags_as_int_r_free_flags(
          ma, test_flag_value="f")
      crys_id = 0
      for l in labels:
        if 'crystal_id' in l:
          crys_id = int(l.split('=')[-1])
          break
      if crys_id > 0 and crystal_id is None:
        label += "%i" %crys_id
      if crystal_id is not None and crys_id > 0 and crys_id != crystal_id:
        continue
      if crys_id not in mtz_crystals:
        mtz_crystals[crys_id] = (
          mtz_object.add_crystal(
            name="crystal_%i" %crys_id,
            project_name="project",
            unit_cell=unit_cell), {})
      crystal, datasets = mtz_crystals[crys_id]
      w_id = 0
      for l in labels:
        if 'wavelength_id' in l:
          w_id = int(l.split('=')[-1])
          break
      if wavelength_id is not None and w_id > 0 and w_id != wavelength_id:
        continue
      if w_id > 0 and wavelength_id is None:
        label += "%i" %w_id
      if w_id not in datasets:
        datasets[w_id] = crystal.add_dataset(
          name="dataset",
          wavelength=0)
      dataset = datasets[w_id]
      if not ma.is_unique_set_under_symmetry():
        if merge_non_unique_under_symmetry:
          print "Warning: merging non-unique data"
          try:
            ma = ma.merge_equivalents().array().customized_copy(
              crystal_symmetry=ma).set_info(ma.info())
          except Sorry, e:
            if ("merge_equivalents_exact: incompatible" in str(e)) :
              raise Sorry(str(e) + " for %s" %ma.info().labels[-1])
            raise
        else:
          n_all = ma.indices().size()
          sel_unique = ma.unique_under_symmetry_selection()
          sel_dup = ~flex.bool(n_all, sel_unique)
          n_duplicate = sel_dup.count(True)
          n_uus = sel_unique.size()
          print "Miller indices not unique under symmetry:", file_name, \
                "(%d redundant indices out of %d)" % (n_all-n_uus, n_all)
          print "Add --merge to command arguments to force merging data."
          return None
          if (show_details_if_error):
            ma.show_comprehensive_summary(prefix="  ")
            ma.map_to_asu().sort().show_array(prefix="  ")
      if(map_to_asu):
        ma = ma.map_to_asu().set_info(ma.info())
      if(remove_systematic_absences):
        ma = ma.remove_systematic_absences()
      ma = ma.select_indices(indices=flex.miller_index(((0,0,0),)),negate=True) \
        .set_info(ma.info()) # Get rid of fake (0,0,0) reflection in some CIFs
      column_labels.add(label)
      dataset.add_miller_array(ma, column_root_label=label)
  return mtz_object

########################################################################
# PHENIX GUI ROUTINES
#
master_phil = iotbx.phil.parse("""
input
  .caption = This program will convert CIF-formatted structure factors (used \
    by the PDB) to an MTZ file.  Other CIF types (restraints, etc.) will be \
    ignored.  Because the data in the PDB often contains mistakes or lacks \
    symmetry, an optional PDB file is strongly recommended.  If you want the \
    program to generate an MTZ file for a specific PDB ID, you may specify \
    that instead of input files, and the CIF and PDB will be fetched from \
    www.rcsb.org.
  .style = caption_img:icons/custom/phenix.reflection_file_editor.png
{
  pdb_id = None
    .type = str
    .short_caption = PDB ID to retrieve
    .input_size = 80
    .style = bold noauto
  cif_file = None
    .type = path
    .short_caption = CIF data file
    .style = bold noauto OnChange:extract_symm_for_cif
  pdb_file = None
    .type = path
    .short_caption = PDB file
    .style = bold noauto file_type:pdb OnChange:extract_symm_for_cif
  wavelength_id = None
    .type = str
    .short_caption = Wavelength ID
    .help = Not required when only one wavelength is present
    .input_size = 120
    .style = noauto
  crystal_id = None
    .type = str
    .short_caption = Crystal ID
    .help = Not required when only one crystal is present
    .input_size = 120
    .style = noauto
}
crystal_symmetry {
  space_group = None
    .type = space_group
    .style = bold
  unit_cell = None
    .type = unit_cell
    .style = bold
}
output_file_name = None
  .type = path
  .style = new_file bold
options {
  use_model = False
    .type = bool
    .short_caption = Use model to help guess data type
    .help = If false, the model will only be used to extract crystal symmetry.
  merge = False
    .type = bool
    .short_caption = Merge non-unique data
  show_details_if_error = True
    .type = bool
    .short_caption = Show data details for some errors
  show_log = True
    .type = bool
}
""")

# TODO replace the old 'run' method
#
# XXX this is still a little unsophisticated with respect to extracting
# crystal symmetry, but it's meant to be run from the Phenix GUI right now.
def run2 (args,
          log=sys.stdout,
          check_params=True,
          params=None) :
  import mmtbx.command_line.fetch_pdb
  parameter_interpreter = master_phil.command_line_argument_interpreter(
    home_scope="")
  pdb_file = None
  cif_file = None
  sources = []
  for arg in args :
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg) :
        pdb_files = os.path.abspath(arg)
      elif arg.endswith(".cif") or arg.endswith(".cif.txt") :
        cif_file = os.path.abspath(arg)
      else :
        try :
          user_phil = iotbx.phil.parse(file_name=arg)
        except RuntimeError :
          print "Unrecognizable file format for %s" % arg
        else :
          sources.append(user_phil)
    else :
      if arg.startswith("--") :
        arg = arg[2:] + "=True"
      try :
        user_phil = parameter_interpreter.process(arg=arg)
        sources.append(user_phil)
      except RuntimeError :
        print "Unrecognizable parameter %s" % arg
  if (params is None) :
    params = master_phil.fetch(sources=sources).extract()
  symm = None
  if (params.input.pdb_id is not None) :
    params.input.pdb_file = mmtbx.command_line.fetch_pdb.run2(
      args=[params.input.pdb_id],
      log=log)
    params.input.cif_file = mmtbx.command_line.fetch_pdb.run2(
      args=["-x", params.input.pdb_id],
      log=log)
    symm = crystal_symmetry_from_any.extract_from(params.input.pdb_file)
    params.crystal_symmetry.space_group = symm.space_group_info()
    params.crystal_symmetry.unit_cell = symm.unit_cell()
    params.input.pdb_id = None
  if check_params :
    validate_params(params)
  if params.output_file_name is None :
    base, ext = os.path.splitext(params.input.cif_file)
    params.output_file_name = os.path.join(os.getcwd(), base + ".mtz")
  if not params.options.use_model :
    params.input.pdb_file = None
  if symm is None :
    assert (type(params.crystal_symmetry.space_group).__name__ ==
            "space_group_info")
    symm = crystal.symmetry(
      space_group_info=params.crystal_symmetry.space_group,
      unit_cell=params.crystal_symmetry.unit_cell)
  n_refl = process_files(
    file_name=params.input.cif_file,
    crystal_symmetry=symm,
    pdb_file_name=params.input.pdb_file,
    output_file_name=params.output_file_name,
    wavelength_id=params.input.wavelength_id,
    crystal_id=params.input.crystal_id,
    show_details_if_error=params.options.show_details_if_error,
    output_r_free_label="FreeR_flag")
  return (params.output_file_name, n_refl)

def validate_params (params) :
  if (params.input.cif_file is None) and (params.input.pdb_id is None) :
    raise Sorry("No CIF file provided!")
  if (params.input.pdb_id is not None) :
    if (params.input.cif_file is not None) :
      raise Sorry("Please specify either a PDB ID or a CIF file, not both.")
    import iotbx.pdb.fetch
    try :
      iotbx.pdb.fetch.validate_pdb_id(params.input.pdb_id)
    except RuntimeError, e :
      raise Sorry(str(e))
  else :
    if ((params.crystal_symmetry.space_group is None) or
        (params.crystal_symmetry.unit_cell is None)) :
      raise Sorry("Crystal symmetry missing or incomplete.")
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.chdir(self.output_dir)
    return run2(args=list(self.args), log=sys.stdout)

def finish_job (results) :
  (mtz_file, n_refl) = results
  if n_refl is None :
    n_refl = 0
  if (mtz_file is not None) and os.path.isfile(mtz_file) :
    return ([("MTZ file", mtz_file)], [("Number of reflections", n_refl)])
  return ([], [])

if(__name__ == "__main__"):
   run(sys.argv[1:])
