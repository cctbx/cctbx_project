"""Convert CIF to MTZ"""
# LIBTBX_SET_DISPATCHER_NAME phenix.cif_as_mtz
# LIBTBX_SET_DISPATCHER_NAME iotbx.cif_as_mtz

from __future__ import absolute_import, division, print_function
from iotbx.option_parser import iotbx_option_parser
from iotbx import crystal_symmetry_from_any
import iotbx.phil
import iotbx.mtz
import iotbx.pdb
from iotbx import cif_mtz_data_labels
from cctbx.array_family import flex
from cctbx import crystal
from libtbx import runtime_utils
from libtbx.utils import Sorry
import libtbx.callbacks # import dependency
import string
import re
import os
import sys
import six

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

def run(args, command_name = "phenix.cif_as_mtz", out=sys.stdout,
       return_as_miller_arrays=False):
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
      .option(None, "--incompatible_flags_to_work_set",
        action="store_true",
        help="When merging place reflections with incompatible flags into the "
             "working set.")
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
      .option("--ignore_bad_sigmas",
          action="store_true",
          help="Set sigmas to None instead of raising an error when bad sigmas "
               "are present.")
      .option("--extend_flags",
          action="store_true",
          help="Extend R-free flags to cover all reflections if necessary.")

    ).process(args=args)
  except Exception as e:
    if(str(e) != "0"): print(str(e))
    sys.exit(0)
  crystal_symmetry = command_line.symmetry
  if(len(command_line.args) > 1):
    print("%d arguments are given from the command line:"% \
      len(command_line.args), command_line.args, file=out)
    raise Sorry("Please specify one reflection cif file.")
  file_name = command_line.args[0]
  if(not os.path.isfile(file_name)):
    raise Sorry("File is not found: %s"%file_name)
  output_r_free_label = command_line.options.output_r_free_label
  if ((not output_r_free_label[0] in string.ascii_uppercase) or
      (re.search(r"[^a-zA-Z0-9_\-]", output_r_free_label))):
    raise Sorry(("%s is not a suitable column label.  MTZ format requires "+
      "an uppercase letter as the first character, and only alphanumeric "+
      "characters or hyphens in the rest of the string.")% output_r_free_label)
  result=process_files(
    file_name=file_name,
    crystal_symmetry=crystal_symmetry,
    output_file_name=command_line.options.output_file_name,
    wavelength_id=command_line.options.wavelength_id,
    crystal_id=command_line.options.crystal_id,
    show_details_if_error=command_line.options.show_details_if_error,
    output_r_free_label=command_line.options.output_r_free_label,
    merge_non_unique_under_symmetry=command_line.options.merge,
    map_to_asu=command_line.options.map_to_asu,
    remove_systematic_absences=command_line.options.remove_systematic_absences,
    incompatible_flags_to_work_set=command_line.options.incompatible_flags_to_work_set,
    return_as_miller_arrays=return_as_miller_arrays,
    ignore_bad_sigmas=command_line.options.ignore_bad_sigmas,
    extend_flags=command_line.options.extend_flags,
    log=out)
  if return_as_miller_arrays:
    return result

def process_files(file_name,
                   crystal_symmetry,
                   output_file_name,
                   wavelength_id,
                   crystal_id,
                   show_details_if_error,
                   output_r_free_label,
                   merge_non_unique_under_symmetry=False,
                   map_to_asu=False,
                   remove_systematic_absences=False,
                   incompatible_flags_to_work_set=False,
                   ignore_bad_sigmas=False,
                   return_as_miller_arrays=False,
                   extend_flags=False,
                   log=sys.stdout):
  mtz_object = extract(
    file_name                       = file_name,
    crystal_symmetry                = crystal_symmetry,
    wavelength_id                   = wavelength_id,
    crystal_id                      = crystal_id,
    show_details_if_error           = show_details_if_error,
    output_r_free_label             = output_r_free_label,
    merge_non_unique_under_symmetry = merge_non_unique_under_symmetry,
    map_to_asu                      = map_to_asu,
    remove_systematic_absences      = remove_systematic_absences,
    incompatible_flags_to_work_set  = incompatible_flags_to_work_set,
    ignore_bad_sigmas               = ignore_bad_sigmas,
    return_as_miller_arrays         = return_as_miller_arrays,
    extend_flags                    = extend_flags,
    log                             = log)

  if return_as_miller_arrays:
    return mtz_object

  if(mtz_object is not None):
    if not output_file_name :
      basename = os.path.basename(file_name)
      if(basename[-4:-3] == "."): output_file_name = basename[:-4]+".mtz"
      elif(basename[-5:-4] == "."): output_file_name = basename[:-5]+".mtz"
      elif(basename.endswith(".ent.gz")): output_file_name=basename[:-7]+".mtz"
      else: output_file_name = basename+".mtz"
    mtz_object.write(file_name = output_file_name)
    return mtz_object.n_reflections()

def get_label(miller_array, output_r_free_label):
  label = None
  for l in miller_array.info().labels:
    if miller_array.anomalous_flag():
      if miller_array.is_xray_amplitude_array():
        label = "F"
      elif miller_array.is_xray_intensity_array():
        label = "I"
      break
    elif ('_meas' in l):
      if miller_array.is_xray_amplitude_array():
        label = "FOBS"
      elif miller_array.is_xray_intensity_array():
        label = "IOBS"
      elif l.endswith(".phase_meas"):
        label = "PHIM"
      break
    elif ("_calc" in l):
      if miller_array.is_xray_amplitude_array():
        label = "FC"
      elif miller_array.is_xray_intensity_array():
        label = "ICALC"
      elif ".F_calc" in l: # cope with _refln.F_calc_au  and _refln.F_calc labels
        label = "FC"
      elif l.endswith(".phase_calc"):
        label = "PHIC"
      break
    elif 'status' in l or '_free' in l:
      label = output_r_free_label
      break
    elif miller_array.is_hendrickson_lattman_array():
      label = "HL"
      break
    elif (miller_array.is_complex_array()):
      if "DELFWT" in l:
        label = "DELFWT"
        break
      elif "FWT" in l:
        label = "FWT"
        break
    elif (miller_array.is_real_array()):
      if (l.endswith( "pdbx_anom_difference")):
        label = "DANO"
        break
      elif (l.endswith(".fom")):
        label = "FOM"
        break
    # as a last resort try find a match in cif_mtz_data_labels dictionary
    label = cif_mtz_data_labels.ccp4_label_from_cif(l)
    if label:
      return label
  return label

def extract(file_name,
            crystal_symmetry,
            wavelength_id,
            crystal_id,
            show_details_if_error,
            output_r_free_label,
            merge_non_unique_under_symmetry,
            map_to_asu,
            remove_systematic_absences,
            all_miller_arrays=None,
            incompatible_flags_to_work_set=False,
            ignore_bad_sigmas=False,
            extend_flags=False,
            return_as_miller_arrays=False,
            log=sys.stdout):
  import iotbx.cif
  from cctbx import miller
  if all_miller_arrays is None:
    base_array_info = miller.array_info(
      crystal_symmetry_from_file=crystal_symmetry)
    all_miller_arrays = iotbx.cif.reader(file_path=file_name).build_miller_arrays(
      base_array_info=base_array_info)
  if (len(all_miller_arrays) == 0):
    raise Sorry("No data arrays were found in this CIF file.  Please make "+
      "sure that the file contains reflection data, rather than the refined "+
      "model.")
  column_labels = set()
  if (extend_flags):
    map_to_asu = True
  # TODO: is all_mille_arrays a dict ? If not change back
  for (data_name, miller_arrays) in six.iteritems(all_miller_arrays):
    for ma in miller_arrays.values():
      other_symmetry = crystal_symmetry
      try:
        crystal_symmetry = other_symmetry.join_symmetry(
          other_symmetry=ma.crystal_symmetry(),
          force=True)
      except AssertionError as e:
        str_e = str(e)
        from six.moves import cStringIO as StringIO
        s = StringIO()
        if "Space group is incompatible with unit cell parameters." in str_e:
          other_symmetry.show_summary(f=s)
          ma.crystal_symmetry().show_summary(f=s)
          str_e += "\n%s" %(s.getvalue())
          raise Sorry(str_e)
        else:
          raise
  if(crystal_symmetry.unit_cell() is None or
     crystal_symmetry.space_group_info() is None):
    raise Sorry(
      "Crystal symmetry is not defined. Please use the --symmetry option.")
  mtz_object = iotbx.mtz.object() \
    .set_title(title="phenix.cif_as_mtz") \
    .set_space_group_info(space_group_info=crystal_symmetry.space_group_info())
  unit_cell=crystal_symmetry.unit_cell()
  mtz_crystals = {}
  mtz_object.set_hkl_base(unit_cell=unit_cell)
  from iotbx.reflection_file_utils import cif_status_flags_as_int_r_free_flags
  # generate list of all reflections (for checking R-free flags)
  from iotbx.reflection_file_utils import make_joined_set
  all_arrays = []
  for (data_name, miller_arrays) in six.iteritems(all_miller_arrays):
    for ma in miller_arrays.values():
      all_arrays.append(ma)
  complete_set = make_joined_set(all_arrays)
  if return_as_miller_arrays:
    miller_array_list=[]
  current_i = -1
  uc = None
  for i, (data_name, miller_arrays) in enumerate(six.iteritems(all_miller_arrays)):
    for ma in miller_arrays.values():
      #ma = ma.customized_copy(
      #  crystal_symmetry=crystal_symmetry).set_info(ma.info())
      if ma._space_group_info is None:
        ma._space_group_info = crystal_symmetry.space_group_info()
      labels = ma.info().labels
      label = get_label(miller_array=ma, output_r_free_label=output_r_free_label)
      if label is None:
        print("Can't determine output label for %s - skipping." % \
          ma.info().label_string(), file=log)
        continue
      elif label.startswith(output_r_free_label):
        ma, _ = cif_status_flags_as_int_r_free_flags(
          ma, test_flag_value="f")
        if isinstance(ma.data(), flex.double):
          data_int = ma.data().iround()
          assert data_int.as_double().all_eq(ma.data())
          ma = ma.customized_copy(data=data_int).set_info(ma.info())
      elif ((ma.is_xray_amplitude_array() or ma.is_xray_intensity_array())
            and isinstance(ma.data(), flex.int)):
        ma = ma.customized_copy(data=ma.data().as_double()).set_info(ma.info())
      crys_id = 0
      for l in labels:
        if 'crystal_id' in l:
          crys_id = int(l.split('=')[-1])
          break
      if crys_id > 0 and crystal_id is None:
        label += "%i" %crys_id
      if crystal_id is not None and crys_id > 0 and crys_id != crystal_id:
        continue

      if ma.unit_cell() is not None: # use symmetry file on the command line if it's None
        unit_cell = ma.unit_cell()

      if crys_id not in mtz_crystals or \
        (i > current_i and unit_cell is not None and uc is not None and unit_cell.parameters() != uc.parameters()):
        # Ensure new mtz crystals are created if miller_array objects have different unit cells
        # Can happen if there are more datasets in the same cif file, like MAD datasets
        uc = unit_cell
        current_i = i
        # Use unique project and crystal names so that MtzGet() in cmtzlib.c picks up individual unit cells
        mtz_crystals[crys_id] = (
          mtz_object.add_crystal(
            name="crystal_%i" %i,
            project_name="project_%i" %i,
            unit_cell =uc), {})
      crystal, datasets = mtz_crystals[crys_id]
      w_id = 0
      for l in labels:
        if 'wavelength_id' in l:
          w_id = int(l.split('=')[-1])
          break
      if wavelength_id is not None and w_id > 0 and w_id != wavelength_id:
        continue
      if w_id > 1 and wavelength_id is None:
        if (label in column_labels):
          label += "%i" %w_id
        #print "label is", label
      if w_id not in datasets:
        wavelength = ma.info().wavelength
        if (wavelength is None):
          wavelength = 0
        datasets[w_id] = crystal.add_dataset(
          name="dataset",
          wavelength=wavelength)
      dataset = datasets[w_id]
      # if all sigmas for an array are set to zero either raise an error, or set sigmas to None
      if ma.sigmas() is not None and (ma.sigmas() == 0).count(False) == 0:
        if ignore_bad_sigmas:
          print("Warning: bad sigmas, setting sigmas to None.", file=log)
          ma.set_sigmas(None)
        else:
          raise Sorry(
  """Bad sigmas: all sigmas are equal to zero.
  Add --ignore_bad_sigmas to command arguments to leave out sigmas from mtz file.""")
      if not ma.is_unique_set_under_symmetry():
        if merge_non_unique_under_symmetry:
          print("Warning: merging non-unique data", file=log)
          if (label.startswith(output_r_free_label)
              and incompatible_flags_to_work_set):
            merging = ma.merge_equivalents(
              incompatible_flags_replacement=0)
            if merging.n_incompatible_flags > 0:
              print("Warning: %i reflections were placed in the working set " \
                    "because of incompatible flags between equivalents." %(
                      merging.n_incompatible_flags), file=log)
          else:
            try:
              merging = ma.merge_equivalents()
            except Sorry as e:
              if ("merge_equivalents_exact: incompatible" in str(e)):
                raise Sorry(str(e) + " for %s" %ma.info().labels[-1] + "\n" +
                  "Add --incompatible_flags_to_work_set to command line "
                  "arguments to place incompatible flags to working set.")
                raise
          ma = merging.array().customized_copy(
            crystal_symmetry=ma).set_info(ma.info())
        elif return_as_miller_arrays: # allow non-unique set
          pass
        else:
          n_all = ma.indices().size()
          sel_unique = ma.unique_under_symmetry_selection()
          sel_dup = ~flex.bool(n_all, sel_unique)
          n_duplicate = sel_dup.count(True)
          n_uus = sel_unique.size()
          msg = (
            "Miller indices not unique under symmetry: " + file_name + \
            "(%d redundant indices out of %d)" % (n_all-n_uus, n_all) +
            "Add --merge to command arguments to force merging data.")
          if (show_details_if_error):
            print(msg, file=log)
            ma.show_comprehensive_summary(prefix="  ")
            ma.map_to_asu().sort().show_array(prefix="  ")
          raise Sorry(msg)
      if(map_to_asu):
        ma = ma.map_to_asu().set_info(ma.info())
      if(remove_systematic_absences):
        ma = ma.remove_systematic_absences()
      if (label.startswith(output_r_free_label) and complete_set is not None):
        n_missing = len(complete_set.lone_set(other=ma).indices())
        if (n_missing > 0):
          if (extend_flags):
            from cctbx import r_free_utils
            # determine flag values
            fvals = list(set(ma.data()))
            print("fvals", fvals, file=log)
            fval = None
            if(len(fvals)==1):
              fval = fvals[0]
            elif(len(fvals)==2):
              f1 = (ma.data()==fvals[0]).count(True)/ma.data().size()
              f2 = (ma.data()==fvals[1]).count(True)/ma.data().size()
              if(f1<f2): fval = fvals[0]
              else:      fval = fvals[1]
            elif(len(fvals)==0):
              fval = None
            else:
              fval = 0
              if(not fval in fvals):
                raise Sorry("Cannot determine free-R flag value.")
            #
            if(fval is not None):
              ma = r_free_utils.extend_flags(
                r_free_flags=ma,
                test_flag_value=fval,
                array_label=label,
                complete_set=complete_set,
                preserve_input_values=True,
                allow_uniform_flags=True,
                log=sys.stdout)
            else:
              ma = None
          else :
            strings = ["%d reflections do not have R-free flags in the "%n_missing,
              "array '%s' - this may "%label,
              "cause problems if you try to use the MTZ file for refinement ",
              "or map calculation.  We recommend that you extend the flags ",
              "to cover all reflections (--extend_flags on the command line)."]
            print("WARNING: ", "\n".join(strings), file=log)
      # Get rid of fake (0,0,0) reflection in some CIFs
      if(ma is not None):
        ma = ma.select_indices(indices=flex.miller_index(((0,0,0),)),
          negate=True).set_info(ma.info())

      if return_as_miller_arrays:
        miller_array_list.append(ma)
        continue  # don't make a dataset

      dec = None
      if ("FWT" in label):
        dec = iotbx.mtz.ccp4_label_decorator()
      column_types = None
      if ("PHI" in label or "PHWT" in label) and (ma.is_real_array()):
        column_types = "P"
      elif (label.startswith("DANO") and ma.is_real_array()):
        if (ma.sigmas() is not None):
          column_types = "DQ"
        else :
          column_types = "D"
      label_base = label
      i = 1
      while label in column_labels:
        label = label_base + "-%i" %(i)
        i += 1
      if(ma is not None):
        column_labels.add(label)
        if("FWT-1" in label): dec=None
        dataset.add_miller_array(ma,
          column_root_label=label,
          label_decorator=dec,
          column_types=column_types)
  if return_as_miller_arrays:
    return miller_array_list
  else:
    return mtz_object

########################################################################
# PHENIX GUI ROUTINES
#
master_phil = iotbx.phil.parse("""
input
  .caption = This program will convert CIF-formatted structure factors (used \
    by the PDB) to an MTZ file.  Other CIF types (restraints, etc.) will be \
    ignored.  Because the data in the PDB often contains mistakes or lacks \
    symmetry, an optional model file is strongly recommended.  If you want the \
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
    .short_caption = Model file
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
include scope libtbx.phil.interface.tracking_params
options {
  merge = False
    .type = bool
    .short_caption = Merge non-unique data
  map_to_asu = True
    .type = bool
    .short_caption = Map HKL indices to ASU
  eliminate_sys_absent = False
    .type = bool
    .short_caption = Remove systematic absences
  show_details_if_error = True
    .type = bool
    .short_caption = Show data details for some errors
  incompatible_flags_to_work_set = False
    .type = bool
    .short_caption = Move incompatible flags to work set
  ignore_bad_sigmas = False
    .type = bool
  extend_flags = False
    .type = bool
    .short_caption = Extend incomplete R-free flags
  show_log = True
    .type = bool
}
""", process_includes=True)

# TODO replace the old 'run' method
#
# XXX this is still a little unsophisticated with respect to extracting
# crystal symmetry, but it's meant to be run from the Phenix GUI right now.
def run2(args,
          log=sys.stdout,
          check_params=True,
          params=None):
  import mmtbx.command_line.fetch_pdb
  libtbx.call_back.set_warning_log(sys.stderr)
  parameter_interpreter = master_phil.command_line_argument_interpreter(
    home_scope="")
  pdb_file = None
  cif_file = None
  sources = []
  for arg in args :
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdb_files = os.path.abspath(arg)
      elif arg.endswith(".cif") or arg.endswith(".cif.txt"):
        cif_file = os.path.abspath(arg)
      else :
        try :
          user_phil = iotbx.phil.parse(file_name=arg)
        except RuntimeError :
          print("Unrecognizable file format for %s" % arg)
        else :
          sources.append(user_phil)
    else :
      if arg.startswith("--"):
        arg = arg[2:] + "=True"
      try :
        user_phil = parameter_interpreter.process(arg=arg)
        sources.append(user_phil)
      except RuntimeError :
        print("Unrecognizable parameter %s" % arg)
  if (params is None):
    params = master_phil.fetch(sources=sources).extract()
  symm = None
  if (params.input.pdb_id is not None):
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
  if symm is None :
    assert (type(params.crystal_symmetry.space_group).__name__ ==
            "space_group_info")
    symm = crystal.symmetry(
      space_group_info=params.crystal_symmetry.space_group,
      unit_cell=params.crystal_symmetry.unit_cell)
  n_refl = process_files(
    file_name=params.input.cif_file,
    crystal_symmetry=symm,
    output_file_name=params.output_file_name,
    wavelength_id=params.input.wavelength_id,
    crystal_id=params.input.crystal_id,
    show_details_if_error=params.options.show_details_if_error,
    output_r_free_label="FreeR_flag",
    merge_non_unique_under_symmetry=params.options.merge,
    map_to_asu=params.options.map_to_asu,
    remove_systematic_absences=params.options.eliminate_sys_absent,
    incompatible_flags_to_work_set=\
      params.options.incompatible_flags_to_work_set,
    ignore_bad_sigmas=params.options.ignore_bad_sigmas,
    extend_flags=params.options.extend_flags)
  return (params.output_file_name, n_refl)

def validate_params(params):
  if params.input.cif_file is None and params.input.pdb_id is None:
    raise Sorry("No CIF file provided!")
  if params.input.cif_file == [] and params.input.pdb_id is None:
    raise Sorry("No structure factors found!")
  if (params.input.pdb_id is not None):
    if (params.input.cif_file is not None):
      raise Sorry("Please specify either a PDB ID or a CIF file, not both.")
    import iotbx.pdb.fetch
    try :
      iotbx.pdb.fetch.validate_pdb_id(params.input.pdb_id)
    except RuntimeError as e :
      raise Sorry(str(e))
  else :
    if ((params.crystal_symmetry.space_group is None) or
        (params.crystal_symmetry.unit_cell is None)):
      raise Sorry("Crystal symmetry missing or incomplete.")
  if (params.output_file_name is not None):
    output_dir = os.path.dirname(params.output_file_name)
    if (not os.path.isdir(output_dir)):
      raise Sorry(("The output directory %s does not exist or is not a "+
        "directory.") % output_dir)
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.chdir(self.output_dir)
    return run2(args=list(self.args), log=sys.stdout)

def finish_job(results):
  (mtz_file, n_refl) = results
  if n_refl is None :
    n_refl = 0
  if (mtz_file is not None) and os.path.isfile(mtz_file):
    return ([("MTZ file", mtz_file)], [("Number of reflections", n_refl)])
  return ([], [])

if(__name__ == "__main__"):
  run(sys.argv[1:])

