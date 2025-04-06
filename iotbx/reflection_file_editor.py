"""GUI tool for manipulating reflection file data
"""
# TODO: confirm old_test_flag_value if ambiguous

from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx.utils import Sorry, null_out, check_if_output_directory_exists
from libtbx import adopt_init_args, slots_getstate_setstate
import warnings
import random
import string
import re
import os
import sys
import six
from six.moves import range
from six.moves import zip

DEBUG = False

# XXX: note that extend=True in the Phenix GUI
master_phil = iotbx.phil.parse("""
show_arrays = False
  .type = bool
  .help = Command-line option, prints out a list of the arrays in each file
  .style = hidden
dry_run = False
  .type = bool
  .help = Print out final configuration and output summary, but don't write \
          the output file
  .style = hidden
verbose = True
  .type = bool
  .help = Print extra debugging information
  .style = hidden
mtz_file
{
  crystal_symmetry
    .caption = By default, the input crystal symmetry will be used for the \
      output file.
    .style = auto_align menu_item
  {
    unit_cell = None
      .type = unit_cell
      .style = bold noauto
    space_group = None
      .type = space_group
      .style = bold noauto
    output_unit_cell = None
      .type = unit_cell
      .style = bold noauto
    output_space_group = None
      .type = space_group
      .style = bold noauto
    change_of_basis = None
      .type = str
      .expert_level = 2
    eliminate_invalid_indices = False
      .type = bool
      .expert_level = 2
    expand_to_p1 = False
      .type = bool
      .short_caption = Expand to P1
      .style = bold
    disable_unit_cell_check = False
      .type = bool
      .style = noauto
      .short_caption = Disable unit cell isomorphism check
    disable_space_group_check = False
      .type = bool
      .style = noauto
    eliminate_sys_absent = True
      .type = bool
      .short_caption = Eliminate systematic absences
  }
  d_max = None
    .type = float
    .short_caption = Low resolution
    .style = bold renderer:draw_hkltools_resolution_widget
  d_min = None
    .type = float
    .short_caption = High resolution
    .style = bold renderer:draw_hkltools_resolution_widget
  wavelength = None
    .type = float
    .short_caption = Wavelength
    .style = bold
  output_file = None
    .type = path
    .short_caption = Output file
    .style = file_type:mtz new_file bold noauto
  include scope libtbx.phil.interface.tracking_params
  resolve_label_conflicts = False
    .type = bool
    .help = Updates label names as necessary to avoid conflicts
    .short_caption = Automatically resolve output label conflicts
    .style = noauto bold
  exclude_reflection = None
    .multiple = True
    .optional = True
    .type = ints(size=3)
    .input_size = 100
    .style = noauto menu_item
  miller_array
    .multiple = True
    .short_caption = Output data array
    .style = fixed auto_align
  {
    file_name = None
      .type = path
      .style = noedit bold
    labels = None
      .type = str
      .short_caption = Array name
      .style = noedit bold
    output_labels = None
      .type = strings
      .optional = True
      .short_caption = Output column labels
      .input_size = 300
      .help = Most Miller arrays have more than one label, and there must be \
              exactly as many new labels as the number of labels in the \
              old array.  Note however that the output labels do not \
              necessarily correspond to the original array name.  (See caveat \
              in Phenix manual about Scalepack files.)
      .style = fixed
    column_root_label = None
      .type = str
      .optional = True
      .short_caption = Base column label
      .help = If specified, this is applied to all columns in the output \
        array with default prefixes and suffixes.  Overrides the \
        output_labels parameter.
      .input_size = 200
    d_min = None
      .type = float
      .short_caption = High resolution
    d_max = None
      .type = float
      .short_caption = Low resolution
    output_as = *auto intensities amplitudes_fw amplitudes
      .type = choice
      .short_caption = Output diffraction data as
      .caption = Automatic Intensities Amplitudes_(run_French-Wilson) \
        Amplitudes_(simple_conversion)
      .help = If the Miller array is amplitudes or intensities, this flag \
        determines the output data type.  If intensities are being converted \
        to amplitudes, this can optionally be done using the French and \
        Wilson treatment to correct weak and negative values.
    anomalous_data = *Auto merged anomalous
      .type = choice
      .caption = automatic_(keep) merge_if_present force_anomalous_output
      .short_caption = Anomalous data
      .help = Optional averaging or generation of Bijvoet mates.  Note that the \
        number of labels expected may change.
    force_type = *auto amplitudes intensities
      .type = choice
      .short_caption = Force observation type
      .help = If the Miller array is amplitudes or intensites, this flag \
        allows the observation type to be set without modifying the data. \
        This is primarily used for structure factors downloaded from the PDB, \
        which sometimes have the data type specified incorrectly.
    scale_max = None
      .type = float
      .short_caption = Scale to maximum value
      .help = Scales data such that the maximum is equal to the given value
    scale_factor = None
      .type = float
      .help = Multiplies data with the given factor
    remove_negatives = False
      .type = bool
      .short_caption = Remove negative values
    massage_intensities = False
      .type = bool
    filter_by_signal_to_noise = None
      .type = float
      .short_caption = Filter by signal-to-noise ratio
    add_b_iso = None
      .type = float
      .short_caption = Add isotropic B-factor
    add_b_aniso = 0 0 0 0 0 0
      .type = floats(size=6)
      .short_caption = Add anisotropic B-factor
    shuffle_values = False
      .type = bool
    reset_values_to = None
      .type = float
      .short_caption = Reset values to
      .help = If defined, all data values will be set to this number.  Sigmas \
        will be left unmodified.  Only applies to single-value floating-point \
        data (I, F, PHI, etc.).
  }
  r_free_flags
    .short_caption = R-free flags generation
    .style = menu_item auto_align box
  {
    generate = True
      .type = bool
      .short_caption = Generate R-free flags if not already present
      .style = bold noauto
    force_generate = False
      .type = bool
      .short_caption = Generate R-free flags even if they are already present
    new_label = FreeR_flag
      .type = str
      .short_caption = Output label for new R-free flags
      .input_size = 160
      .style = bold noauto
    include scope cctbx.r_free_utils.generate_r_free_params_str
    random_seed = None
      .type = int
      .short_caption = Seed for random number generator
      .expert_level = 2
    extend = None
      .type = bool
      .short_caption = Extend existing R-free array(s) to full resolution range
      .style = bold noauto
    old_test_flag_value = None
      .type = int
      .short_caption = Original test flag value
      .help = Overrides automatic guess of test flag value from existing set. \
        This will usually be 1 for reflection files generated by Phenix, and \
        0 for reflection files from CCP4.  Do not change unless you're sure \
        you know what flag to use!
      .expert_level = 2
    export_for_ccp4 = False
      .type = bool
      .short_caption = Convert R-free flags to CCP4 convention
      .help = If True, R-free flags expressed as boolean values (or 0 and 1) \
        will be converted to random integers from 0 to 19, where 0 denotes \
        the test set.  Phenix will work with either convention, but most CCP4 \
        programs expect the test set to be 0.
      .expert_level = 2
      .style = noauto
    preserve_input_values = True
      .type = bool
      .short_caption = Preserve original flag values
      .help = If True, CCP4-style R-free flags will be left as random \
        integers instead of being converted to a boolean array.  This option \
        is not compatible with the option to export flags to the CCP4 \
        convention.
      .style = noauto
    warn_if_all_same_value = True
      .type = bool
      .short_caption = Warn if R-free flags are all the same value
    adjust_fraction = False
      .type = bool
      .short_caption = Adjust test set size to specified fraction
      .help = If True, the R-free flags will be resized as necessary.  This \
        may be useful when the current set is too large or too small, but \
        needs to be preserved (at least in part).  The target fraction will \
        be relative to all possible reflections, even those not measured in \
        the input.  Note that this option is not compatible with preserving \
        input values.
    d_eps = 0.0001
      .type = float
      .short_caption = Resolution buffer
      .expert_level = 2
    relative_to_complete_set = False
      .type = bool
      .short_caption = Generate R-free flags relative to complete set
    remediate_mismatches = False
      .type = bool
      .short_caption = Remediate mismatching Friedel/Bijvoet pairs
      .help = If True, flags that differ between F+ and F- will be moved into \
        the test set for both reflections.  This option is incompatible with \
        preserving the input values as integers.
  }
}""", process_includes=True)

class array_input(slots_getstate_setstate):
  __slots__ = ["miller_arrays", "file_names", "array_types"]
  def __init__(self, params, input_files=None):
    from iotbx import file_reader
    if (input_files is None):
      input_files = {}
    self.miller_arrays = []
    self.file_names = []
    self.array_types = []
    for i_array, array_params in enumerate(params.mtz_file.miller_array):
      if (array_params.file_name is None):
        raise Sorry("Missing file name for array %d (labels=%s)" %
          (i_array+1, str(array_params.labels)))
      elif (not os.path.isfile(array_params.file_name)):
        raise Sorry("The path '%s' does not exist or is not a file." %
          array_params.file_name)
      input_file = input_files.get(array_params.file_name)
      if input_file is None :
        input_file = file_reader.any_file(array_params.file_name,
          force_type="hkl")
      is_mtz = (input_file.file_object.file_type() == "ccp4_mtz")
      if input_file is None :
        input_file = file_reader.any_file(array_params.file_name)
        input_file.assert_file_type("hkl")
        input_files[input_file.file_name] = input_file
      found = False
      for miller_array in input_file.file_object.as_miller_arrays():
        array_info = miller_array.info()
        label_string = array_info.label_string()
        if label_string == array_params.labels :
          self.miller_arrays.append(miller_array)

          self.file_names.append(input_file.file_name)
          found = True
          if is_mtz :
            array_type = get_original_array_types(input_file,
              original_labels=array_info.labels)
            self.array_types.append(array_type)
          else :
            self.array_types.append(None)
      if not found :
        raise Sorry("Couldn't find the Miller array %s in file %s!" %
                    (array_params.labels, array_params.file_name))

  def resolve_unit_cell(self):
    unit_cells = []
    for any_array_type in [False, True] :
      for array in self.miller_arrays :
        if array.is_experimental_data():
          cs = array.crystal_symmetry()
          if (cs is not None):
            uc = cs.unit_cell()
            unit_cells.append(uc)
      if (len(unit_cells) > 0):
        break
    if (len(unit_cells) == 0):
      raise Sorry("No unit cell found in inputs - please specify explicitly.")
    return unit_cells

class process_arrays(object):
  def __init__(self, params, input_files=None, inputs_object=None,
      log=sys.stderr, accumulation_callback=None, symmetry_callback=None):
    if (input_files is None):
      input_files = {}
    adopt_init_args(self, locals())
    validate_params(params)
    r_free_params = params.mtz_file.r_free_flags
    import cctbx.miller
    from cctbx import r_free_utils
    from cctbx import crystal
    from scitbx.array_family import flex

    #-------------------------------------------------------------------
    # COLLECT ARRAYS
    if (inputs_object is None):
      inputs_object = array_input(params=params, input_files=input_files)
    miller_arrays = inputs_object.miller_arrays
    file_names = inputs_object.file_names
    array_types = inputs_object.array_types
    if params.show_arrays :
      shown_files = []
      for file_name, miller_array in zip(file_names, miller_arrays):
        if not file_name in shown_files :
          print("%s:" % file_name, file=log)
          shown_files.append(file_name)
        print("  %s" % miller_array.info().label_string(), file=log)
      return

    labels = ["H", "K", "L"]
    label_files = [None, None, None]
    self.created_r_free = False
    have_r_free_array = False
    self.final_arrays = []
    self.mtz_dataset = None
    self.wavelength = params.mtz_file.wavelength
    if (self.wavelength is None):
      file_wavelengths = {}
      for miller_array in miller_arrays :
        if (not miller_array.is_xray_data_array()) : continue
        info = miller_array.info()
        if (info is not None):
          if (not info.source in file_wavelengths):
            if (info.wavelength is not None):
              try :
                file_wavelengths[info.source] = float(info.wavelength)
              except ValueError as e :
                print("Warning: bad wavelength '%s'" % info.wavelength, file=log)
      all_wavelengths = set([ w for f, w in six.iteritems(file_wavelengths) ])
      if (len(all_wavelengths) == 1):
        self.wavelength = all_wavelengths.pop()
      elif (len(all_wavelengths) > 1):
        filtered_wavelengths = set([])
        # MTZ files typically have the wavelength set to 1.0 if it was not
        # previously specified, so we ignore these.
        for file_name, wavelength in six.iteritems(file_wavelengths):
          if (file_name.endswith(".mtz") and wavelength in [0.0, 1.0]):
            continue
          filtered_wavelengths.add(wavelength)
        if (len(filtered_wavelengths) == 1):
          self.wavelength = filtered_wavelengths.pop()
        elif (len(filtered_wavelengths) > 1):
          raise Sorry(("Multiple wavelengths present in input experimental "+
            "data arrays: %s.  Please specify the wavelength parameter "+
            "explicitly.") %
            ", ".join([ "%g"%x for x in sorted(list(filtered_wavelengths)) ]))

    #-------------------------------------------------------------------
    # SYMMETRY SETUP
    change_symmetry = change_point_group = False
    input_space_group = params.mtz_file.crystal_symmetry.space_group
    output_space_group = params.mtz_file.crystal_symmetry.output_space_group
    if (output_space_group is not None):
      output_sg = params.mtz_file.crystal_symmetry.output_space_group.group()
      input_point_group = input_space_group.group().build_derived_point_group()
      output_point_group = \
        output_space_group.group().build_derived_point_group()
      pg_number_in = input_point_group.type().number()
      pg_number_out = output_point_group.type().number()
      change_symmetry = True
      if (pg_number_out != pg_number_in):
        change_point_group = True
        print("Will expand to P1 symmetry before merging.", file=log)
    else :
      output_sg = params.mtz_file.crystal_symmetry.space_group.group()
    if params.mtz_file.crystal_symmetry.output_unit_cell is not None :
      output_uc = params.mtz_file.crystal_symmetry.output_unit_cell
      change_symmetry = True
    else :
      output_uc = params.mtz_file.crystal_symmetry.unit_cell
    if params.mtz_file.crystal_symmetry.expand_to_p1 and change_symmetry :
      raise Sorry("Output unit cell and space group must be undefined if "+
          "expand_to_p1 is True.")
    input_symm = crystal.symmetry(
      unit_cell=params.mtz_file.crystal_symmetry.unit_cell,
      space_group_info=params.mtz_file.crystal_symmetry.space_group,
      assert_is_compatible_unit_cell=False,
      force_compatible_unit_cell=False)
    if (not input_symm.is_compatible_unit_cell()):
      raise Sorry(("Input unit cell %s is incompatible with the specified "+
        "space group (%s).") % (str(params.mtz_file.crystal_symmetry.unit_cell),
          str(params.mtz_file.crystal_symmetry.space_group)))
    derived_sg = input_symm.space_group().build_derived_point_group()
    output_symm = crystal.symmetry(
      unit_cell=output_uc,
      space_group=output_sg,
      assert_is_compatible_unit_cell=False,
      force_compatible_unit_cell=False)
    if (not output_symm.is_compatible_unit_cell()):
      raise Sorry(("Output unit cell %s is incompatible with the specified "+
        "space group (%s).") % (str(output_uc), str(output_sg)))

    # Resolution limits
    (d_max, d_min) = get_best_resolution(miller_arrays, input_symm)
    if (d_max is None) and (params.mtz_file.d_max is None):
      raise Sorry("No low-resolution cutoff could be found in the "+
        "parameters or input file(s); you need to explicitly set this value "+
        "for the program to run.")
    if d_min is None :
      if params.mtz_file.d_min is None :
        raise Sorry("No high-resolution cutoff could be found in the "+
          "parameters or input file(s); you need to explicitly set this "+
          "value for the program to run.")
    # XXX resolution limits are used even if outside the range used by the
    # input data, in case we want to generate extra R-free flags for future
    # use
    if (params.mtz_file.d_max is not None):
      d_max = params.mtz_file.d_max
    if (params.mtz_file.d_min is not None):
      d_min = params.mtz_file.d_min
    if r_free_params.random_seed is not None :
      random.seed(r_free_params.random_seed)
      flex.set_random_seed(r_free_params.random_seed)

    #-------------------------------------------------------------------
    # MAIN LOOP
    i = 0
    r_free_arrays = []
    for (array_params, file_name, miller_array) in \
        zip(params.mtz_file.miller_array, file_names, miller_arrays):
      info = miller_array.info()
      array_name = "%s:%s" % (file_name, array_params.labels)

      #-----------------------------------------------------------------
      # APPLY SYMMETRY
      array_sg = miller_array.space_group()
      array_uc = miller_array.unit_cell()
      ignore_sg = params.mtz_file.crystal_symmetry.disable_space_group_check
      ignore_uc = params.mtz_file.crystal_symmetry.disable_unit_cell_check
      if (array_sg is not None) and (not ignore_sg):
        if array_sg.build_derived_point_group() != derived_sg :
          raise Sorry(("The point group for the Miller array %s (%s) does "+
            "not match the point group of the overall space group (%s).") %
            (array_name, str(array_sg), str(input_symm.space_group())))
      if (array_uc is not None) and (not ignore_uc):
        if not array_uc.is_similar_to(input_symm.unit_cell()):
          raise Sorry(("The unit cell for the Miller array %s (%s) is "+
            "significantly different than the output unit cell (%s).  You "+
            "can ignore this by setting disable_unit_cell_check=True (in "+
            "the GUI, enable \"Disable unit cell isomorphism check\").") %
            (array_name, str(array_uc), str(input_symm.unit_cell())))
      miller_array = miller_array.customized_copy(
        crystal_symmetry=input_symm).map_to_asu()
      if params.mtz_file.crystal_symmetry.expand_to_p1 :
        miller_array = miller_array.expand_to_p1()
      elif change_symmetry :
        sg_number = miller_array.space_group_info().type().number()
        if (change_point_group) and (sg_number != 1):
          miller_array = miller_array.expand_to_p1()
        miller_array = miller_array.customized_copy(
          crystal_symmetry=output_symm)
      if not miller_array.is_unique_set_under_symmetry():
        if miller_array.is_integer_array() and not is_rfree_array(miller_array,info):
          raise Sorry(("The data in %s cannot be merged because they are in "+
            "integer format.  If you wish to change symmetry (or the input "+
            "data are unmerged), you must omit this array.  (Note also that "+
            "merging will fail for R-free flags if the flags for symmetry-"+
            "related reflections are not identical.)") % array_name)
        miller_array = miller_array.merge_equivalents().array()

      if DEBUG :
        print("  Adjusted size:  %d" % miller_array.data().size(), file=log)
        if miller_array.sigmas() is not None :
          print("         sigmas:  %d" % miller_array.sigmas().size(), file=log)

      #-----------------------------------------------------------------
      # CHANGE OF BASIS
      # this will actually be done at the last minute - we just check for
      # incompatible options here
      if params.mtz_file.crystal_symmetry.change_of_basis is not None :
        if change_symmetry :
          raise Sorry("You may not change symmetry when change_of_basis is "+
            "defined.")

      output_array = modify_array(
        miller_array=miller_array,
        array_params=array_params,
        array_name=array_name,
        array_info=info,
        log=log,
        verbose=params.verbose)
      if([params.mtz_file.d_min,params.mtz_file.d_max].count(None)==0):
        if(params.mtz_file.d_min > params.mtz_file.d_max):
          msg="High resolution cutoff %s is larger than low resolution %s."
          raise Sorry(msg%(
            str(params.mtz_file.d_min), str(params.mtz_file.d_max)))
      output_array = output_array.resolution_filter(
        d_min=params.mtz_file.d_min,
        d_max=params.mtz_file.d_max)
      if (len(params.mtz_file.exclude_reflection) > 0):
        for hkl in params.mtz_file.exclude_reflection :
          output_array = output_array.delete_index(hkl)

      #-----------------------------------------------------------------
      # OUTPUT
      assert isinstance(output_array, cctbx.miller.array)
      default_label = array_params.column_root_label
      output_labels = array_params.output_labels
      if (default_label is not None):
        output_labels = None
      else :
        if i <= 25:
          default_label = 2 * string.ascii_uppercase[i]
        else:
          i1 = int(i//25)
          i2 = i - i1 * 25
          if i2 > 25:
            raise Sorry("Maximum of 625 columns in reflection file editor")

          default_label = 2 * string.ascii_uppercase[i1]  + \
             2* string.ascii_uppercase[i2]
      column_types = None
      import iotbx.mtz
      default_types = iotbx.mtz.default_column_types(output_array)
      if (array_types[i] is not None):
        if output_array.is_xray_amplitude_array():
          array_types[i] = re.sub("J", "F", array_types[i])
        elif output_array.is_xray_intensity_array():
          array_types[i] = re.sub("F", "J", array_types[i])
        if len(default_types) == len(array_types[i]):
          print("Recovering original column types %s" % array_types[i], file=log)
          column_types = array_types[i]
      if (output_array.data().size() == 0):
        raise Sorry("The array %s:%s ended up empty.  Please check the "+
          "resolution cutoffs to make sure they do not exclude all data "+
          "from the input file.  If you think the parameters were correct, "+
          "this is probably a bug; please contact help@phenix-online.org "+
          "with a description of the problem.")
      if is_rfree_array(output_array, info):
        r_free_arrays.append((output_array, info, output_labels,
          array_params.column_root_label, file_name))
      else :
        if DEBUG :
          print("  Final size:    %d" % output_array.data().size(), file=log)
          if output_array.sigmas() is not None :
            print("      sigmas:    %d" % output_array.sigmas().size(), file=log)
        try:
          self.add_array_to_mtz_dataset(
            output_array=output_array,
            default_label=default_label,
            column_types=column_types,
            out=log)
          ok_array = True
        except Exception as e:
          ok_array = False
          print("Skipping array '%s' which cannot be converted to MTZ" %(
            output_array.info), file = log)
        if ok_array:
          if (output_labels is not None):
            for label in output_labels :
              labels.append(label)
              label_files.append(file_name)
          else :
            n_cols = len(default_types)
            if (output_array.anomalous_flag() and
                not output_array.is_xray_reconstructed_amplitude_array()):
              n_cols *= 2
            for k in range(n_cols):
              labels.append(None)
              label_files.append(file_name)
          self.final_arrays.append(output_array)

      i += 1

    #-------------------------------------------------------------------
    # EXISTING R-FREE ARRAYS
    if len(r_free_arrays) > 0 :
      from iotbx.reflection_file_utils import get_r_free_flags_scores, \
        make_joined_set
      have_r_free_array = True
      combined_set = complete_set = None
      if len(self.final_arrays) > 0 :
        eps = r_free_params.d_eps
        combined_set = make_joined_set(self.final_arrays)
        complete_set = combined_set.complete_set(d_min=d_min-eps,
          d_max=d_max+eps)
        # XXX big hack
        missing = combined_set.lone_set(other=complete_set)
        if (missing.size() > 0):
          complete_set = complete_set.concatenate(other=missing)
      if (len(self.final_arrays) > 1):
        warnings.warn("Multiple Miller arrays are already present in this "+
          "file; the R-free flags will be generated based on the total "+
          "of reflections in all arrays combined.  If you want the fraction "+
          "of test set reflections to be relative to a specific array, you "+
          "should run the editor with that array separately first.",
          UserWarning)
      if (r_free_params.relative_to_complete_set):
        combined_set = complete_set
      elif (combined_set is not None):
        check_and_warn_about_incomplete_r_free_flags(combined_set)
      i = 0
      for (new_array, info, output_labels, root_label, file_name) in \
          r_free_arrays :
        # XXX this is important for guessing the right flag when dealing
        # with CCP4-style files, primarily when the flag values are not
        # very evenly distributed
        new_array.set_info(info)
        test_flag_value = None
        flag_scores = get_r_free_flags_scores(miller_arrays=[new_array],
          test_flag_value=r_free_params.old_test_flag_value)
        test_flag_value = flag_scores.test_flag_values[0]
        if (test_flag_value is None):
          if (r_free_params.old_test_flag_value is not None):
            test_flag_value = r_free_params.old_test_flag_value
          elif ((r_free_params.warn_if_all_same_value) or
                (r_free_params.extend)):
            raise Sorry(("The data in %s:%s appear to be R-free flags, but "+
              "a suitable test flag value (usually 1 or 0) could not be "+
              "automatically determined.  This may indicate that the flags "+
              "are uniform, which is not suitable for refinement; it can "+
              "also happen when there are exactly 3 different values used. "+
              " If this is not the case, you may specify the test flag "+
              "value manually by clicking the button labeled \"R-free flags "+
              "generation...\" and entering the value to use under "+
              "\"Original test flag value\".  Alternately, unchecking the box "+
              "\"Warn if R-free flags are all the same value\" will skip this "+
              "step, but you will not be able to extend the flags to higher "+
              "resolution.") % (file_name, info.label_string()))
        if (r_free_params.preserve_input_values):
          assert (not r_free_params.remediate_mismatches)
          r_free_flags = new_array
        else :
          new_data = (new_array.data()==test_flag_value)
          assert isinstance(new_data, flex.bool)
          r_free_flags = new_array.array(data=new_data)
        r_free_flags = r_free_flags.map_to_asu()
        generate_bijvoet_mates = (array_params.anomalous_data=="anomalous")
        if not r_free_flags.is_unique_set_under_symmetry():
          r_free_flags = r_free_flags.merge_equivalents().array()
        if (r_free_flags.anomalous_flag()):
          if (r_free_params.remediate_mismatches):
            print("Remediating any mismatched flags for Friedel mates...", file=log)
            r_free_flags = r_free_utils.remediate_mismatches(
              array=r_free_flags,
              log=log)
          r_free_flags = r_free_flags.average_bijvoet_mates()
          if (output_labels is not None) and (len(output_labels) != 1):
            assert (not combined_set.anomalous_flag())
            # XXX can't do this operation on a miller set - will expand the
            # r-free flags later
            generate_bijvoet_mates = True
        if (r_free_params.adjust_fraction):
          print("Resizing test set in %s" % array_name, file=log)
          r_free_as_bool = get_r_free_as_bool(r_free_flags,
            test_flag_value)
          r_free_flags = r_free_utils.adjust_fraction(
            miller_array=r_free_as_bool,
            fraction=r_free_params.fraction,
            log=log)
        if (r_free_params.extend):
          r_free_flags = r_free_utils.extend_flags(
            r_free_flags=r_free_flags,
            test_flag_value=test_flag_value,
            array_label=array_name,
            complete_set=combined_set,
            accumulation_callback=accumulation_callback,
            preserve_input_values=r_free_params.preserve_input_values,
            d_max=d_max,
            d_min=d_min,
            log=log)
        output_array = r_free_flags
        if (generate_bijvoet_mates):
          output_array = output_array.generate_bijvoet_mates()
        if (len(params.mtz_file.exclude_reflection) > 0):
          for hkl in params.mtz_file.exclude_reflection :
            output_array = output_array.delete_index(hkl)
        if (r_free_params.export_for_ccp4):
          print("%s: converting to CCP4 convention" % array_name, file=log)
          output_array = export_r_free_flags(
            miller_array=output_array,
            test_flag_value=True)
        default_label = root_label
        if (default_label is None):
          default_label = "A" + string.ascii_uppercase[i+1]
        self.add_array_to_mtz_dataset(
          output_array=output_array,
          default_label=default_label,
          column_types="I",
          out=log)
        if (output_labels is not None):
          validate_output_labels(
            miller_array=output_array,
            array_params=array_params,
            array_name=array_name,
            output_labels=output_labels)
          for label in output_labels :
            labels.append(label)
            label_files.append(file_name)
        else :
          if (output_array.anomalous_flag()):
            labels.extend([None,None])
            label_files.extend([file_name,file_name])
          else :
            labels.append(None)
            label_files.append(file_name)
        self.final_arrays.append(output_array)
        i += 1

    #-------------------------------------------------------------------
    # NEW R-FREE ARRAY
    if ((r_free_params.generate and not have_r_free_array) or
        r_free_params.force_generate):
      if (len(self.final_arrays) > 1):
        warnings.warn("Multiple Miller arrays are already present in this "+
          "file; the R-free flags will be generated based on the total "+
          "of reflections in all arrays combined.  If you want the fraction "+
          "of test set reflections to be relative to a specific array, you "+
          "should run the editor with that array separately first.",
          UserWarning)
      from iotbx.reflection_file_utils import make_joined_set
      combined_set = make_joined_set(self.final_arrays)
      complete_set = combined_set.complete_set(
        d_min=d_min-r_free_params.d_eps,
        d_max=d_max+r_free_params.d_eps)
      if (r_free_params.relative_to_complete_set):
        combined_set = complete_set
      else :
        check_and_warn_about_incomplete_r_free_flags(combined_set)
      # XXX this used to generate the test set from the actually complete set,
      # but the fraction free needs to be relative to the
      new_r_free_array = combined_set.generate_r_free_flags(
        fraction=r_free_params.fraction,
        max_free=r_free_params.max_free,
        lattice_symmetry_max_delta=r_free_params.lattice_symmetry_max_delta,
        use_lattice_symmetry=r_free_params.use_lattice_symmetry,
        use_dataman_shells=r_free_params.use_dataman_shells,
        n_shells=r_free_params.n_shells)
      if r_free_params.new_label is None or r_free_params.new_label == "" :
        r_free_params.new_label = "FreeR_flag"
      if r_free_params.export_for_ccp4 :
        print("%s: converting to CCP4 convention" % array_name, file=log)
        output_array = export_r_free_flags(
          miller_array=new_r_free_array,
          test_flag_value=True)
      else:
        output_array = new_r_free_array
      if (len(params.mtz_file.exclude_reflection) > 0):
        for hkl in params.mtz_file.exclude_reflection :
          output_array = output_array.delete_index(hkl)
      self.add_array_to_mtz_dataset(
        output_array=output_array,
        default_label="ZZZZ",
        column_types="I",
        out=log)
      labels.append(r_free_params.new_label)
      label_files.append("(new array)")
      self.created_r_free = True

    #-------------------------------------------------------------------
    # RE-LABEL COLUMNS
    mtz_object = self.mtz_dataset.mtz_object()
    self.label_changes = []
    self.mtz_object = mtz_object
    if not len(labels) == mtz_object.n_columns():
      print("\n".join([ "LABEL: %s" % label for label in labels ]), file=log)
      self.show(out=log)
      raise Sorry("The number of output labels does not match the final "+
        "number of labels in the MTZ file.  Details have been printed to the "+
        "console.")
    i = 0
    used = dict([ (label, 0) for label in labels ])
    invalid_chars = re.compile(r"[^A-Za-z0-9_\-+\(\)]")
    for column in self.mtz_object.columns():
      if (labels[i] is not None) and (column.label() != labels[i]):
        label = labels[i]
        original_label = label
        if invalid_chars.search(label) is not None :
          raise Sorry(("Invalid label '%s'.  Output labels may only contain "+
            "alphanumeric characters, underscore, plus and minus signs, or "+
            "parentheses.")
            % label)
        if used[label] > 0 :
          if params.mtz_file.resolve_label_conflicts :
            if label.endswith("(+)") or label.endswith("(-)"):
              label = label[0:-3] + ("_%d" % (used[label]+1)) + label[-3:]
            else :
              label += "_%d" % (used[labels[i]] + 1)
            if used.get(label,None) is not None:  # this was already there
              found = False
              for k in range(1,10000):
                if label.endswith("(+)") or label.endswith("(-)"):
                  new_label = \
                     label[0:-3] + ("_%d" % (used[label]+k)) + label[-3:]
                else :
                  new_label = label + "_%d" % (used[labels[i]] + k)
                if used.get(new_label,None) is None: # got it
                  label = new_label
                  used[label] = 0
                  found = True
                  break
              if not found:
                raise Sorry("Unable to resolve multiple similar labels")
            self.label_changes.append((label_files[i], labels[i], label))
          else :
            raise Sorry(("Duplicate column label '%s'.  Specify "+
              "resolve_label_conflicts=True to automatically generate "+
              "non-redundant labels, or edit the parameters to provide your"+
              "own choice of labels.") % labels[i])
        try :
          column.set_label(label)
        except RuntimeError as e :
          if ("new_label is used already" in str(e)):
            col_names = [ col.label() for col in mtz_object.columns() ]
            raise RuntimeError(("Duplicate column label '%s': current labels "+
              "are %s; user-specified output labels are %s.") %
              (label, " ".join(col_names), " ".join(labels)))
        else :
          used[original_label] += 1
      i += 1

  def add_array_to_mtz_dataset(self, output_array, default_label,
      column_types, out=sys.stdout):
    # apply change of basis here
    if (self.params.mtz_file.crystal_symmetry.change_of_basis is not None):
      output_array, cb_op = output_array.apply_change_of_basis(
        change_of_basis=self.params.mtz_file.crystal_symmetry.change_of_basis,
        eliminate_invalid_indices=\
          self.params.mtz_file.crystal_symmetry.eliminate_invalid_indices,
        out=out)
    if (self.params.mtz_file.crystal_symmetry.eliminate_sys_absent):
      output_array = output_array.eliminate_sys_absent()
    if self.mtz_dataset is None :
      self.mtz_dataset = output_array.as_mtz_dataset(
        column_root_label=default_label,
        column_types=column_types,
        wavelength=self.wavelength)
    else :
      self.mtz_dataset.add_miller_array(
        miller_array=output_array,
        column_root_label=default_label,
        column_types=column_types)

  def show(self, out=sys.stdout):
    if self.mtz_object is not None :
      print("", file=out)
      print(("=" * 20) + " Summary of output file " + ("=" * 20), file=out)
      self.mtz_object.show_summary(out=out, prefix="  ")
      print("", file=out)

  def finish(self):
    assert self.mtz_object is not None
    if self.params.verbose :
      self.show(out=self.log)
    self.mtz_object.write(file_name=self.params.mtz_file.output_file)
    n_refl = self.mtz_object.n_reflections()
    del self.mtz_object
    self.mtz_object = None
    return n_refl

#-----------------------------------------------------------------------
# TODO get rid of these two (need to make sure they aren't imported elsewhere)
def get_r_free_stats(*args, **kwds):
  from cctbx import r_free_utils
  return r_free_utils.get_r_free_stats(*args, **kwds)

def get_r_free_as_bool(*args, **kwds):
  from cctbx import r_free_utils
  return r_free_utils.get_r_free_as_bool(*args, **kwds)

def get_best_resolution(miller_arrays, input_symm=None):
  best_d_min = None
  best_d_max = None
  for array in miller_arrays :
    array_symm = array.crystal_symmetry()
    if ((None in [array_symm.space_group(), array_symm.unit_cell()]) and
        (input_symm is not None)):
      array = array.customized_copy(crystal_symmetry=input_symm)
    try :
      (d_max, d_min) = array.d_max_min()
      if best_d_max is None or d_max > best_d_max :
        best_d_max = d_max
      if best_d_min is None or d_min < best_d_min :
        best_d_min = d_min
    except Exception as e :
      pass
  return (best_d_max, best_d_min)

def is_rfree_array(miller_array, array_info):
  from iotbx import reflection_file_utils
  return ((miller_array.is_integer_array() or
           miller_array.is_bool_array()) and
          reflection_file_utils.looks_like_r_free_flags_info(array_info))

def export_r_free_flags(miller_array, test_flag_value):
  from cctbx import r_free_utils
  new_flags = r_free_utils.export_r_free_flags_for_ccp4(
    flags=miller_array.data(),
    test_flag_value=test_flag_value)
  return miller_array.customized_copy(data=new_flags)

def get_original_array_types(input_file, original_labels):
  array_types = ""
  mtz_file = input_file.file_object.file_content()
  mtz_columns = mtz_file.column_labels()
  mtz_types = mtz_file.column_types()
  mtz_crossref = dict(zip(mtz_columns, mtz_types))
  for label in original_labels :
    array_types += mtz_crossref[label]
  return array_types

def guess_array_output_labels(miller_array):
  info = miller_array.info()
  assert info is not None
  labels = info.labels
  output_labels = labels
  if (labels in [["i_obs","sigma"], ["Intensity+-","SigmaI+-"]]):
    if miller_array.anomalous_flag():
      output_labels = ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    else :
      output_labels = ["I", "SIGI"]
  elif (miller_array.is_xray_reconstructed_amplitude_array()):
    output_labels = ["F", "SIGF", "DANO", "SIGDANO", "ISYM"]
  elif ((miller_array.is_xray_amplitude_array() or
         miller_array.is_xray_intensity_array()) and
        miller_array.anomalous_flag()):
    if (len(labels) == 2) and (miller_array.sigmas() is not None):
      output_labels = [ "%s(+)" % labels[0],  # data(+)
                        "%s(+)" % labels[1],  # sigma(+)
                        "%s(-)" % labels[0],  # data(-)
                        "%s(-)" % labels[1] ] # sigma(-)
    elif (len(labels) == 1) and (miller_array.sigmas() is None):
      output_labels = [ "%s(+)" % labels[0], "%s(-)" % labels[0] ]
    elif (miller_array.anomalous_flag()) and (len(labels) == 5): # mmCIF
      output_labels = ["I(+)","SIGI(+)","I(-)","SIGI(-)"]
    elif (not miller_array.anomalous_flag()) and (len(labels) == 3): # mmCIF
      output_labels = ["I","SIGI"]

  # Catch general mmCIF inputs and remove leading label (it is dataset name)
  #  and remove anything before "." because those are not allowed in output
  n_expected = get_number_of_expected_columns(miller_array,
    raise_sorry_on_errors = False)
  if output_labels and (len(output_labels) == n_expected + 1):
    output_labels = edit_mmcif_output_labels(output_labels,
      n_expected = n_expected)
  return output_labels

def edit_mmcif_output_labels(output_labels, n_expected = None):
  assert n_expected and len(output_labels) == n_expected + 1
  output_labels = output_labels[1:]
  new_output_labels = []
  for o in output_labels:
    if o.find(".") > -1:
      o = o.split(".")[-1]
    new_output_labels.append(o)
  return new_output_labels

def modify_array(
    miller_array,
    array_name,
    array_params,
    array_info=None,
    verbose=True,
    debug=False,
    log=sys.stdout):
  """
  Perform various manipulations on a Miller array.  This can be applied to
  any data type, although certain options are limited to experimental data.
  """
  from scitbx.array_family import flex
  output_labels = array_params.output_labels
  if (output_labels is None) and (array_params.column_root_label is None):
    raise Sorry("Missing output labels for %s!" % array_name)
  if (array_params.column_root_label is not None):
    output_labels = None
  labels_base = re.sub(",merged$", "", array_params.labels)
  input_labels = labels_base.split(",")
  if not None in [array_params.scale_factor, array_params.scale_max] :
    raise Sorry("The parameters scale_factor and scale_max are " +
      "mutually exclusive.")
  if not False in [array_params.remove_negatives,
                   array_params.massage_intensities] :
    raise Sorry("The parameters remove_negatives and massage_intensities "+
      "are mutually exclusive.")
  if debug :
    print("  Starting size:  %d" % miller_array.data().size(), file=log)
    if miller_array.sigmas() is not None :
      print("         sigmas:  %d" % miller_array.sigmas().size(), file=log)
  if verbose :
    print("Processing %s" % array_name, file=log)
  if array_params.d_max is not None and array_params.d_max <= 0 :
    array_params.d_max = None
  if array_params.d_min is not None and array_params.d_min <= 0 :
    array_params.d_min = None
  output_labels = array_params.output_labels
  if (output_labels is None) and (array_params.column_root_label is None):
    raise Sorry("Missing output labels for %s!" % array_name)
  if (array_params.column_root_label is not None):
    output_labels = None
  wavelength = getattr(array_info, "wavelength", None)
  if not None in [array_params.scale_factor, array_params.scale_max] :
    raise Sorry("The parameters scale_factor and scale_max are " +
      "mutually exclusive.")
  if not False in [array_params.remove_negatives,
                   array_params.massage_intensities] :
    raise Sorry("The parameters remove_negatives and massage_intensities "+
      "are mutually exclusive.")
  if debug :
    print("  Starting size:  %d" % miller_array.data().size(), file=log)
    if miller_array.sigmas() is not None :
      print("         sigmas:  %d" % miller_array.sigmas().size(), file=log)
  if debug :
    print("  Resolution before resolution filter: %.2f - %.2f" % (
      miller_array.d_max_min()), file=log)
  # go ahead and perform the array-specific cutoff
  miller_array = miller_array.resolution_filter(
    d_min=array_params.d_min,
    d_max=array_params.d_max)
  if debug :
    print("              after resolution filter: %.2f - %.2f" % (
      miller_array.d_max_min()), file=log)
    print("  Truncated size: %d" % miller_array.data().size(), file=log)
    if miller_array.sigmas() is not None :
      frint >> log, "          sigmas: %d" % miller_array.sigmas().size()
  # anomalous manipulation
  if (miller_array.anomalous_flag() and
      array_params.anomalous_data == "merged"):
    print(("Converting array %s from anomalous to non-anomalous." %
                   array_name), file=log)
    if (not miller_array.is_xray_intensity_array()):
      miller_array = miller_array.average_bijvoet_mates()
      if miller_array.is_xray_reconstructed_amplitude_array():
        miller_array.set_observation_type_xray_amplitude()
    else :
      miller_array = miller_array.average_bijvoet_mates()
      miller_array.set_observation_type_xray_intensity()
  elif ((not miller_array.anomalous_flag()) and
        (array_params.anomalous_data == "anomalous")):
    print("Generating Bijvoet mates for %s" % array_name, file=log)
    miller_array = miller_array.generate_bijvoet_mates()
  # scale factors
  if (array_params.scale_max is not None):
    print(("Scaling %s such that the maximum value is: %.6g" %
                   (array_name, array_params.scale_max)), file=log)
    miller_array = miller_array.apply_scaling(target_max=array_params.scale_max)
  elif (array_params.scale_factor is not None):
    print(("Multiplying data in %s with the factor: %.6g" %
                   (array_name, array_params.scale_factor)), file=log)
    miller_array = miller_array.apply_scaling(factor=array_params.scale_factor)
  # Since this function has many built-in consistency checks, we always run
  # it even for non-experimental arrays
  miller_array = modify_experimental_data_array(
    miller_array=miller_array,
    array_name=array_name,
    array_params=array_params,
    log=log)
  # Information removal for running control experiments - normally these
  # would be used on experimental data, but there is no reason why they can't
  # be applied to other array types.
  if (array_params.shuffle_values):
    print("Shuffling values for %s" % array_name, file=log)
    combined_array = None
    tmp_array = miller_array.deep_copy()
    tmp_array.setup_binner(n_bins=min(100, tmp_array.indices().size()//10))
    for i_bin in tmp_array.binner().range_used():
      bin_sel = tmp_array.binner().selection(i_bin)
      bin_array = tmp_array.select(bin_sel)
      perm = flex.random_permutation(bin_array.data().size())
      sigmas = bin_array.sigmas()
      if (sigmas is not None):
        sigmas = sigmas.select(perm)
      data = bin_array.data().select(perm)
      bin_array = bin_array.customized_copy(data=data, sigmas=sigmas)
      if (combined_array is None):
        combined_array = bin_array.deep_copy()
      else :
        combined_array = combined_array.concatenate(bin_array)
    if (combined_array.indices().size() != miller_array.indices().size()):
      raise RuntimeError("Array size changed: %d versus %d" %
        (combined_array.indices().size(), miller_array.indices().size()))
    miller_array = combined_array
  if (array_params.reset_values_to):
    if (not miller_array.is_real_array() and
        not miller_array.is_integer_array()):
      raise Sorry("Resetting the values for %s is not permitted." %
        array_name)
    print("Resetting values for %s to %g" % (array_name,
      array_params.reset_values_to), file=log)
    data = miller_array.data().deep_copy()
    new_value = array_params.reset_values_to
    if miller_array.is_integer_array():
      new_value = int(new_value)
    data.fill(new_value)
    miller_array = miller_array.customized_copy(data=data)
  if (array_params.force_type != "auto"):
    if (not miller_array.is_xray_data_array()):
      raise Sorry(("You may only override the output observation type for "+
        "amplitudes or intensities - the data in %s are unsupported.") %
        array_name)
    if (array_params.force_type == "amplitudes"):
      miller_array = miller_array.set_observation_type_xray_amplitude()
    elif (array_params.force_type == "intensities"):
      miller_array = miller_array.set_observation_type_xray_intensity()
  if (not is_rfree_array(miller_array, array_info)):
    if (array_params.column_root_label is not None):
      validate_column_root_label(
        miller_array=miller_array,
        array_name=array_name,
        root_label=array_params.column_root_label)
    elif (output_labels is not None):
      validate_output_labels(
        miller_array=miller_array,
        array_params=array_params,
        array_name=array_name,
        output_labels=output_labels)
  return miller_array

def modify_experimental_data_array(
    miller_array,
    array_name,
    array_params,
    log=sys.stdout):
  """
  Manipulations common to amplitude and intensity arrays only.
  """
  # negative intensity/amplitude remediation
  if (array_params.remove_negatives):
    if (miller_array.is_real_array()):
      print("Removing negatives from %s" % array_name, file=log)
      miller_array = miller_array.select(miller_array.data() > 0)
      if (miller_array.sigmas() is not None):
        miller_array = miller_array.select(miller_array.sigmas() > 0)
    else :
      raise Sorry("remove_negatives not applicable to %s." % array_name)
  elif (array_params.massage_intensities):
    if (miller_array.is_xray_intensity_array()):
      if array_params.output_as == "amplitudes" :
        miller_array = miller_array.enforce_positive_amplitudes()
      else :
        raise Sorry(("You must output %s as amplitudes to use the "+
          "massage_intensities option.") % array_name)
    else :
      raise Sorry("The parameter massage_intensities is only valid for "+
        "X-ray intensity arrays.")
  # I/sigma filtering
  if (array_params.filter_by_signal_to_noise is not None):
    if (not miller_array.is_xray_data_array()):
      raise Sorry(("Filtering by signal-to-noise is only supported for "+
        "amplitudes or intensities (failed on %s).") % array_name)
    elif (array_params.filter_by_signal_to_noise <= 0):
      raise Sorry(("A value greater than zero is required for the "+
        "cutoff for filtering by signal to noise ratio (failed array: %s).") %
        array_name)
    sigmas = miller_array.sigmas()
    if (sigmas is None):
      raise Sorry(("Sigma values must be defined to filter by signal "+
        "to noise ratio (failed on %s).") % array_name)
    elif (not sigmas.all_ne(0.0)):
      # XXX should it just remove these too?
      raise Sorry(("The sigma values for the array %s include one or "+
        "more zeros - filtering by signal to noise not supported.") %
        array_name)
    data = miller_array.data()
    miller_array = miller_array.select(
      (data / sigmas) > array_params.filter_by_signal_to_noise)
  # apply B-factors
  # leave the default as [0]*6 to make the format clear, but reset to
  # None if unchanged
  if (array_params.add_b_aniso == [0,0,0,0,0,0]):
    array_params.add_b_aniso = None
  if ((array_params.add_b_iso is not None) or
      (array_params.add_b_aniso is not None)):
    if (not miller_array.is_real_array() and
        not miller_array.is_complex_array()):
      raise Sorry(("Applying a B-factor to the data in %s is not "+
        "permitted.") % array_name)
    if (array_params.add_b_iso is not None):
      miller_array = miller_array.apply_debye_waller_factors(
        b_iso=array_params.add_b_iso,
        apply_to_sigmas=True)
    if (array_params.add_b_aniso is not None):
      miller_array = miller_array.apply_debye_waller_factors(
        b_cart=array_params.add_b_aniso,
        apply_to_sigmas=True)
  # data type manipulation
  if miller_array.is_xray_intensity_array():
    if (array_params.output_as in ["amplitudes", "amplitudes_fw"]):
      if (array_params.output_as == "amplitudes"):
        miller_array = miller_array.f_sq_as_f()
      else :
        from cctbx import french_wilson
        miller_array = french_wilson.french_wilson_scale(
          miller_array=miller_array,
          log=log)
      miller_array.set_observation_type_xray_amplitude()
  elif miller_array.is_xray_amplitude_array():
    if array_params.output_as == "intensities" :
      miller_array = miller_array.f_as_f_sq()
      miller_array.set_observation_type_xray_intensity()
  return miller_array

def validate_output_labels(
    miller_array,
    array_name,
    array_params,
    output_labels):
  """
  Check output labels for consistency with selected options (done before the
  array undergoes most modifications).
  """
  def raise_sorry_if_wrong_number_of_labels(n_expected):
    assert (n_expected is not None)
    msg_extra = "If you are "+ \
        "unsure which labels to provide, use the simpler column root label "+ \
        "instead, which will be expanded as necessary. (Note that the " + \
        "number of input labels from non-MTZ file formats does not " + \
        "necessarily correspond to the expected number of output labels.)"
    n_used = len(output_labels)
    msg_expected = "%d are" % n_expected
    if (n_expected == 1):
      msg_expected = "%d is" % n_expected
    msg_used = "labels \"%s\" have" % " ".join(output_labels)
    if (n_used == 1):
      msg_used = "label \"%s\" has" % " ".join(output_labels)
    if (n_expected > n_used):
      raise Sorry(("You have not specified enough MTZ column labels for the "+
        "array %s; only the %s been specified, but %s "+
        "required for the output array after processing.  " + msg_extra) %
        (array_name, msg_used, msg_expected))
    elif (n_expected < n_used):
      raise Sorry(("You have specified too many column labels for the array "+
        "%s; the %s been specified, but only %s are allowed "+
        "for the output array after processing.  " + msg_extra) %
        (array_name, msg_used, msg_expected))
  n_expected = get_number_of_expected_columns(miller_array,
    array_params = array_params,
    array_name = array_name,
    output_labels = output_labels,
    raise_sorry_on_errors = True)
  raise_sorry_if_wrong_number_of_labels(n_expected)
  if miller_array.is_xray_data_array():
    validate_column_root_label(
      miller_array=miller_array,
      root_label = output_labels[0].upper(),
      array_name=array_name)
  return True

def get_number_of_expected_columns(miller_array,
    output_labels = None,
    array_params = None,
    array_name = None,
    raise_sorry_on_errors = True):
  if raise_sorry_on_errors:
    assert (output_labels is not None) and (array_params is not None) and (
     array_name is not None)
  n_expected = None
  if (miller_array.is_xray_reconstructed_amplitude_array()):
    # FIXME this needs to be handled better - but it should at least
    # catch files from CCP4 data processing
    if (raise_sorry_on_errors) and ((len(output_labels) != 5) and
        (not array_params.anomalous_data == "merged")):
      raise Sorry(("The array %s will be output as "+
        "amplitudes and anomalous differences with sigmas, plus ISYM. "+
        "Five columns will be written, but %d labels were specified.") %
        (array_name, len(output_labels)))
    n_expected = 5
  elif (miller_array.anomalous_flag()):
    if miller_array.is_real_array():
      if (miller_array.sigmas() is not None):
        n_expected = 4
      else :
        n_expected = 2
    elif miller_array.is_complex_array():
      assert miller_array.sigmas() is None
      n_expected = 4
    elif miller_array.is_hendrickson_lattman_array():
      n_expected = 8
    else :
      n_expected = 2
  else :
    if miller_array.is_real_array():
      if (miller_array.sigmas() is not None):
        n_expected = 2
      else :
        n_expected = 1
    elif miller_array.is_complex_array():
      if (raise_sorry_on_errors) and (
           miller_array.sigmas() is not None):
        raise RuntimeError("Combination of sigmas and complex data not allowed for array %s" % array_name)
      n_expected = 2
    elif miller_array.is_hendrickson_lattman_array():
      n_expected = 4
    else :
      n_expected = 1
  return n_expected


def validate_column_root_label(miller_array, root_label, array_name):
  """
  Check the root MTZ label for a Miller array to ensure consistency with
  data type - this is done after the array has been processed.
  """
  root_label = root_label.upper()
  if (miller_array.is_xray_intensity_array()):
    if (not root_label.startswith("I")):
      raise Sorry(("The column label prefix for the array '%s' "+
        "is inconsistent with the output array type (intensities). "+
        "Please use 'I' (either case) as the first character in the "+
        "label.") % (array_name))
  elif (miller_array.is_xray_amplitude_array()):
    if (not root_label.startswith("F")):
      raise Sorry(("The specified column label prefix for the array '%s' "+
        "is inconsistent with the output array type (amplitudes). "+
        "Please use 'F' (either case) as the first character in the "+
        "label.") % (array_name))
  else :
    if (root_label == "I") or (root_label == "F"):
      raise Sorry(("You have specified the column label prefix '%s' "+
        "for the array '%s', which is neither intensities nor "+
        "amplitudes; the base labels 'I' and 'F' are reserved "+
        "for these array data types .") %
        (root_label, array_name))
  return True

# XXX the requirement for defined crystal symmetry in phil input files is
# problematic for automation, where labels and operations are very
# standardized.  so I added this to identify unique symmetry information
# in the collected input files.
def collect_symmetries(file_names):
  from iotbx import crystal_symmetry_from_any
  file_names = set(file_names)
  file_symmetries = []
  for file_name in file_names :
    symm = crystal_symmetry_from_any.extract_from(file_name)
    if (symm is not None):
      file_symmetries.append(symm)
  return file_symmetries

def resolve_symmetry(file_symmetries, current_space_group, current_unit_cell):
  space_groups = []
  unit_cells = []
  for symm in file_symmetries :
    if (symm is not None):
      space_group = symm.space_group_info()
      unit_cell = symm.unit_cell()
      if (space_group is not None) and (unit_cell is not None):
        if (current_space_group is None):
          group = space_group.group()
          for other_sg in space_groups :
            other_group = other_sg.group()
            if (other_sg.type().number() != group.type().number()):
              raise Sorry("Ambiguous space group information in input files - "+
                "please specify symmetry parameters.")
        if (current_unit_cell is None):
          for other_uc in unit_cells :
            if (not other_uc.is_similar_to(unit_cell, 0.001, 0.1)):
              raise Sorry("Ambiguous unit cell information in input files - "+
                "please specify symmetry parameters.")
        space_groups.append(space_group)
        unit_cells.append(unit_cell)
  consensus_space_group = current_space_group
  consensus_unit_cell = current_unit_cell
  if (len(space_groups) > 0) and (current_space_group is None):
    consensus_space_group = space_groups[0]
  if (len(unit_cells) > 0) and (current_unit_cell is None):
    consensus_unit_cell = unit_cells[0]
  return (consensus_space_group, consensus_unit_cell)

def check_and_warn_about_incomplete_r_free_flags(combined_set):
  completeness = combined_set.completeness()
  if (completeness < 0.99):
    warnings.warn(("The arrays in the input file are incomplete "+
      "(%.1f%% of reflections missing), so the newly generated R-free "+
      "flags will be incomplete as well.  This will not cause a problem "+
      "for refinement, but may result in inconsistent flags later on if "+
      "you collect additional data.  You can disable this warning by "+
      "clicking the box labeled \"Generate flags relative to maximum "+
      "complete set\" in the R-free options dialog window.") %
      (100 * (1-completeness)), UserWarning)

def usage(out=sys.stdout, attributes_level=0):
  print("""
# usage: iotbx.reflection_file_editor [file1.mtz ...] [parameters.eff]
#            --help      (print this message)
#            --details   (show parameter help strings)
# Dumping default parameters:
""", file=out)
  master_phil.show(out=out, attributes_level=attributes_level)

def generate_params(file_name, miller_array, include_resolution=False):
  param_str = """mtz_file.miller_array {
  file_name = %s
  labels = %s
""" % (file_name, miller_array.info().label_string())
  output_labels = guess_array_output_labels(miller_array)
  param_str += "  output_labels = " + " ".join(output_labels)
  if include_resolution :
    try :
      (d_max, d_min) = miller_array.d_max_min()
      # FIXME gross gross gross!
      d_max += 0.001
      d_min -= 0.001
      param_str += """  d_max = %.5f\n  d_min = %.5f\n""" % (d_max, d_min)
    except Exception :
      pass
  param_str += "}"
  return param_str

def validate_params(params):
  if (len(params.mtz_file.miller_array) == 0):
    raise Sorry("No Miller arrays have been selected for the output file.")
  elif len(params.mtz_file.miller_array) > 625 :
    raise Sorry("Only 625 or fewer arrays may be used.")
  if None in [params.mtz_file.crystal_symmetry.space_group,
              params.mtz_file.crystal_symmetry.unit_cell] :
    raise Sorry("Missing or incomplete symmetry information.")
  if (params.mtz_file.r_free_flags.preserve_input_values):
    if (not (0 < params.mtz_file.r_free_flags.fraction < 0.5)):
      raise Sorry("The R-free flags fraction must be greater than zero and "+
        "less than 0.5.")
    if (params.mtz_file.r_free_flags.export_for_ccp4):
      raise Sorry("r_free_flags.preserve_input_values and "+
        "r_free_flags.export_for_ccp4 may not be used together.")
    if (params.mtz_file.r_free_flags.adjust_fraction):
      raise Sorry("Preserving input values of R-free flags is not supported "+
        "when resizing a test set to the specified fraction.")
    if (params.mtz_file.r_free_flags.remediate_mismatches):
      raise Sorry("Preserving input values of R-free flags is not supported "+
        "when correcting mismatched Friedel/Bijvoet pairs.")
  check_if_output_directory_exists(file_name=params.mtz_file.output_file)

#-----------------------------------------------------------------------
def run(args, out=sys.stdout):
  from iotbx import file_reader
  crystal_symmetry_from_pdb = None
  crystal_symmetries_from_hkl = []
  reflection_files = []
  reflection_file_names = []
  user_phil = []
  all_arrays = []
  if len(args) == 0 :
    usage()
  interpreter = master_phil.command_line_argument_interpreter(
    home_scope=None)
  for arg in args :
    if arg in ["--help", "--options", "--details"] :
      usage(attributes_level=args.count("--details"))
      return True
    elif arg in ["-q", "--quiet"] :
      out = null_out()
    elif os.path.isfile(arg):
      full_path = os.path.abspath(arg)
      try :
        file_phil = iotbx.phil.parse(file_name=full_path)
      except Exception :
        pass
      else :
        user_phil.append(file_phil)
        continue
      input_file = file_reader.any_file(full_path)
      if input_file.file_type == "hkl" :
        miller_arrays = input_file.file_server.miller_arrays
        for array in miller_arrays :
          symm = array.crystal_symmetry()
          if symm is not None :
            crystal_symmetries_from_hkl.append(symm)
            break
        reflection_files.append(input_file)
        reflection_file_names.append(os.path.abspath(input_file.file_name))
      elif input_file.file_type == "pdb" :
        symm = input_file.file_object.crystal_symmetry()
        if symm is not None :
          crystal_symmetry_from_pdb = symm
    else :
      if arg.startswith("--"):
        arg = arg[2:] + "=True"
      try :
        cmdline_phil = interpreter.process(arg=arg)
      except Exception :
        pass
      else :
        user_phil.append(cmdline_phil)
  input_files = {}
  for input_file in reflection_files :
    input_files[input_file.file_name] = input_file
    file_arrays = input_file.file_object.as_miller_arrays()
    params = []
    for miller_array in file_arrays :
      params.append(generate_params(input_file.file_name, miller_array))
    file_params_str = "\n".join(params)
    user_phil.append(iotbx.phil.parse(file_params_str))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()
  if crystal_symmetry_from_pdb is not None :
    params.mtz_file.crystal_symmetry.space_group = \
      crystal_symmetry_from_pdb.space_group_info()
    params.mtz_file.crystal_symmetry.unit_cell = \
      crystal_symmetry_from_pdb.unit_cell()
  elif None in [params.mtz_file.crystal_symmetry.space_group,
                params.mtz_file.crystal_symmetry.unit_cell] :
    # extract symmetry information from all reflection files
    all_hkl_file_names = [ p.file_name for p in params.mtz_file.miller_array ]
    need_symmetry_for_files = []
    for file_name in all_hkl_file_names :
      file_name = os.path.abspath(file_name)
      if (not file_name in reflection_file_names):
        need_symmetry_for_files.append(file_name)
    crystal_symmetries_from_hkl.extend(
      collect_symmetries(need_symmetry_for_files))
    (space_group, unit_cell) = resolve_symmetry(
      file_symmetries=crystal_symmetries_from_hkl,
      current_space_group=params.mtz_file.crystal_symmetry.space_group,
      current_unit_cell=params.mtz_file.crystal_symmetry.unit_cell)
    params.mtz_file.crystal_symmetry.space_group = space_group
    params.mtz_file.crystal_symmetry.unit_cell = unit_cell
  if params.mtz_file.output_file is None :
    n = 0
    for file_name in os.listdir(os.getcwd()):
      if file_name.startswith("reflections_") and file_name.endswith(".mtz"):
        n += 1
    params.mtz_file.output_file = "reflections_%d.mtz" % (n+1)
  params.mtz_file.output_file = os.path.abspath(params.mtz_file.output_file)
  process = process_arrays(params, input_files=input_files, log=out)
  if params.dry_run :
    print("# showing final parameters", file=out)
    master_phil.format(python_object=params).show(out=out)
    if params.verbose :
      process.show(out=out)
    return process
  process.finish()
  print("Data written to %s" % params.mtz_file.output_file, file=out)
  return process

if __name__ == "__main__" :
  run(sys.argv[1:])

#---end
