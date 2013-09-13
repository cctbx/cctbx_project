
# TODO: confirm old_test_flag_value if ambiguous

from __future__ import division
import iotbx.phil
from libtbx.utils import Sorry, null_out, check_if_output_directory_exists
from libtbx import adopt_init_args
import warnings
import random
import string
import re
import os
import sys

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
  }
}""", process_includes=True)

class process_arrays (object) :
  def __init__ (self, params, input_files=None, log=sys.stderr,
      accumulation_callback=None, symmetry_callback=None) :
    if (input_files is None) :
      input_files = {}
    adopt_init_args(self, locals())
    validate_params(params)
    r_free_params = params.mtz_file.r_free_flags
    from iotbx import file_reader
    import cctbx.miller
    from cctbx import r_free_utils
    from cctbx import crystal
    from scitbx.array_family import flex

    #-------------------------------------------------------------------
    # COLLECT ARRAYS
    miller_arrays = []
    file_names = []
    array_types = []
    for i_array, array_params in enumerate(params.mtz_file.miller_array) :
      if (array_params.file_name is None) :
        raise Sorry("Missing file name for array %d (labels=%s)" %
          (i_array+1, str(array_params.labels)))
      elif (not os.path.isfile(array_params.file_name)) :
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
      for miller_array in input_file.file_object.as_miller_arrays() :
        array_info = miller_array.info()
        label_string = array_info.label_string()
        if label_string == array_params.labels :
          miller_arrays.append(miller_array)
          file_names.append(input_file.file_name)
          found = True
          if is_mtz :
            array_type = get_original_array_types(input_file,
              original_labels=array_info.labels)
            array_types.append(array_type)
          else :
            array_types.append(None)
      if not found :
        raise Sorry("Couldn't fine the Miller array %s in file %s!" %
                    (array_params.labels, array_params.file_name))
    if params.show_arrays :
      shown_files = []
      for file_name, miller_array in zip(file_names, miller_arrays) :
        if not file_name in shown_files :
          print >> log, "%s:" % file_name
          shown_files.append(file_name)
        print >> log, "  %s" % miller_array.info().label_string()
      return

    labels = ["H", "K", "L"]
    label_files = [None, None, None]
    self.created_r_free = False
    have_r_free_array = False
    self.final_arrays = []
    self.mtz_dataset = None

    #-------------------------------------------------------------------
    # SYMMETRY SETUP
    change_symmetry = False
    if params.mtz_file.crystal_symmetry.output_space_group is not None :
      output_sg = params.mtz_file.crystal_symmetry.output_space_group.group()
      change_symmetry = True
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
    if (not input_symm.is_compatible_unit_cell()) :
      raise Sorry(("Input unit cell %s is incompatible with the specified "+
        "space group (%s).") % (str(params.mtz_file.crystal_symmetry.unit_cell),
          str(params.mtz_file.crystal_symmetry.space_group)))
    derived_sg = input_symm.space_group().build_derived_point_group()
    output_symm = crystal.symmetry(
      unit_cell=output_uc,
      space_group=output_sg,
      assert_is_compatible_unit_cell=False,
      force_compatible_unit_cell=False)
    if (not output_symm.is_compatible_unit_cell()) :
      raise Sorry(("Output unit cell %s is incompatible with the specified "+
        "space group (%s).") % (str(output_uc), str(output_sg)))

    # Resolution limits
    (d_max, d_min) = get_best_resolution(miller_arrays, input_symm)
    if (d_max is None) and (params.mtz_file.d_max is None) :
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
    if (params.mtz_file.d_max is not None) :
      d_max = params.mtz_file.d_max
    if (params.mtz_file.d_min is not None) :
      d_min = params.mtz_file.d_min
    if r_free_params.random_seed is not None :
      random.seed(r_free_params.random_seed)
      flex.set_random_seed(r_free_params.random_seed)

    #-------------------------------------------------------------------
    # MAIN LOOP
    i = 0
    special_labels = ["i_obs,sigma", "Intensity+-,SigmaI+-"]
    r_free_arrays = []
    for (array_params, file_name, miller_array) in \
        zip(params.mtz_file.miller_array, file_names, miller_arrays) :
      array_name = "%s:%s" % (file_name, array_params.labels)
      if params.verbose :
        print >> log, "Processing %s" % array_name
      if array_params.d_max is not None and array_params.d_max <= 0 :
        array_params.d_max = None
      if array_params.d_min is not None and array_params.d_min <= 0 :
        array_params.d_min = None
      output_array = None # this will eventually be the final processed array
      output_labels = array_params.output_labels
      if (output_labels is None) and (array_params.column_root_label is None) :
        raise Sorry("Missing output labels for %s!" % array_name)
      if (array_params.column_root_label is not None) :
        output_labels = None
      info = miller_array.info()
      if not None in [array_params.scale_factor, array_params.scale_max] :
        raise Sorry("The parameters scale_factor and scale_max are " +
          "mutually exclusive.")
      if not False in [array_params.remove_negatives,
                       array_params.massage_intensities] :
        raise Sorry("The parameters remove_negatives and massage_intensities "+
          "are mutually exclusive.")
      if DEBUG :
        print >> log, "  Starting size:  %d" % miller_array.data().size()
        if miller_array.sigmas() is not None :
          print >> log, "         sigmas:  %d" % miller_array.sigmas().size()

      is_experimental_data = (miller_array.is_xray_amplitude_array() or
        miller_array.is_xray_intensity_array() or
        miller_array.is_xray_reconstructed_amplitude_array())

      #-----------------------------------------------------------------
      # OUTPUT LABELS SANITY CHECK
      labels_base = re.sub(",merged$", "", array_params.labels)
      input_labels = labels_base.split(",")
      if (output_labels is not None) :
        if (miller_array.anomalous_flag()) :
          if (array_params.anomalous_data == "merged") :
            if (labels_base in special_labels) :
              if (len(output_labels) != 2) :
                raise Sorry(("There are too many output labels for the array "+
                  "%s, which is being converted to non-anomalous data. "+
                  "Labels such as I,SIGI are appropriate, or F,SIGF if you "+
                  "are converting the array to amplitudes.  (Current output "+
                  "labels: %s)") % (array_name, ",".join(output_labels)))
            elif (len(output_labels) == len(input_labels)) :
              raise Sorry(("There are too many output labels for the array "
                "%s, which is being converted to non-anomalous data. "+
                "The total number of columns will be halved in the output "+
                "array, and the labels should not have trailing (+) or (-).") %
                array_name)
          elif (labels_base in special_labels) and (len(output_labels) != 4) :
            raise Sorry(("There are not enough output labels for the array "+
              "%s. For Scalepack or d*TREK files containing anomalous "+
              "data, you must specify exactly four column labels (e.g. "+
              "I(+),SIGI(+),I(-),SIGI(-), or F(+),SIGF(+),F(-),SIGF(-) "+
              "if you are converting the array to amplitudes).") %
              array_name)
        elif (array_params.anomalous_data == "anomalous") :
          if (len(output_labels) != (len(input_labels) * 2)) :
            raise Sorry(("There are not enough output labels for the array "
              "%s, which is being converted to anomalous data. "+
              "The total number of columns will be doubled in the output "+
              "array, and the labels should be repeated with '(+)' and "+
              "'(-)'.") % array_name)
        if (miller_array.is_xray_reconstructed_amplitude_array()) :
          # FIXME this needs to be handled better - but it should at least
          # catch files from CCP4 data processing
          if ((len(output_labels) != 5) and
              (not array_params.anomalous_data == "merged")) :
            raise Sorry(("The array in %s with labels %s will be output as "+
              "amplitudes and anomalous differences with sigmas, plus ISYM. "+
              "Five columns will be written, but %d labels were specified.") %
              (file_name, array_params.labels, len(output_labels)))

      #-----------------------------------------------------------------
      # APPLY SYMMETRY
      array_sg = miller_array.space_group()
      array_uc = miller_array.unit_cell()
      ignore_sg = params.mtz_file.crystal_symmetry.disable_space_group_check
      ignore_uc = params.mtz_file.crystal_symmetry.disable_unit_cell_check
      if (array_sg is not None) and (not ignore_sg) :
        if array_sg.build_derived_point_group() != derived_sg :
          raise Sorry(("The point group for the Miller array %s (%s) does "+
            "not match the point group of the overall space group (%s).") %
            (array_name, str(array_sg), str(input_symm.space_group())))
      if (array_uc is not None) and (not ignore_uc) :
        if not array_uc.is_similar_to(input_symm.unit_cell()) :
          raise Sorry(("The unit cell for the Miller array %s (%s) is "+
            "significantly different than the output unit cell (%s).  You "+
            "can ignore this by setting disable_unit_cell_check=True (in "+
            "the GUI, enable \"Disable unit cell isomorphism check\").") %
            (array_name, str(array_uc), str(input_symm.unit_cell())))
      new_array = miller_array.customized_copy(
        crystal_symmetry=input_symm).map_to_asu()
      if params.mtz_file.crystal_symmetry.expand_to_p1 :
        new_array = new_array.expand_to_p1()
      elif change_symmetry :
        new_array = new_array.expand_to_p1()
        new_array = new_array.customized_copy(crystal_symmetry=output_symm)
      if not new_array.is_unique_set_under_symmetry() :
        if new_array.is_integer_array() and not is_rfree_array(new_array,info):
          raise Sorry(("The data in %s cannot be merged because they are in "+
            "integer format.  If you wish to change symmetry (or the input "+
            "data are unmerged), you must omit this array.  (Note also that "+
            "merging will fail for R-free flags if the flags for symmetry-"+
            "related reflections are not identical.)") % array_name)
        new_array = new_array.merge_equivalents().array()

      if DEBUG :
        print >> log, "  Adjusted size:  %d" % new_array.data().size()
        if miller_array.sigmas() is not None :
          print >> log, "         sigmas:  %d" % new_array.sigmas().size()

      #-----------------------------------------------------------------
      # CHANGE OF BASIS
      # this will actually be done at the last minute - we just check for
      # incompatible options here
      if params.mtz_file.crystal_symmetry.change_of_basis is not None :
        if change_symmetry :
          raise Sorry("You may not change symmetry when change_of_basis is "+
            "defined.")

      #-----------------------------------------------------------------
      # OTHER FILTERING
      if DEBUG :
        print >> log, "  Resolution before resolution filter: %.2f - %.2f" % (
          new_array.d_max_min())
      # first the array-specific cutoff
      new_array = new_array.resolution_filter(
        d_min=array_params.d_min,
        d_max=array_params.d_max)
      if DEBUG :
        print >> log, "              after resolution filter: %.2f - %.2f" % (
          new_array.d_max_min())
      # now apply the global cutoff
      new_array = new_array.resolution_filter(
          d_min=params.mtz_file.d_min,
          d_max=params.mtz_file.d_max)
      if DEBUG :
        print >> log, "  Truncated size: %d" % new_array.data().size()
        if new_array.sigmas() is not None :
          print >> log, "          sigmas: %d" % new_array.sigmas().size()
      if (new_array.anomalous_flag() and
          array_params.anomalous_data == "merged") :
        print >> log, ("Converting array %s from anomalous to non-anomalous." %
                       array_name)
        if (not new_array.is_xray_intensity_array()) :
          new_array = new_array.average_bijvoet_mates()
        else :
          new_array = new_array.f_sq_as_f()
          new_array = new_array.average_bijvoet_mates()
          new_array = new_array.f_as_f_sq()
          new_array.set_observation_type_xray_intensity()
      elif ((not new_array.anomalous_flag()) and
            (array_params.anomalous_data == "anomalous")) :
        new_array = new_array.generate_bijvoet_mates()
      if (array_params.scale_max is not None) :
        print >> log, ("Scaling %s such that the maximum value is: %.6g" %
                       (array_name, array_params.scale_max))
        new_array = new_array.apply_scaling(target_max=array_params.scale_max)
      elif (array_params.scale_factor is not None) :
        print >> log, ("Multiplying data in %s with the factor: %.6g" %
                       (array_name, array_params.scale_factor))
        new_array = new_array.apply_scaling(factor=array_params.scale_factor)
      if (array_params.remove_negatives) :
        if (new_array.is_real_array()) :
          print >> log, "Removing negatives from %s" % array_name
          new_array = new_array.select(new_array.data() > 0)
          if (new_array.sigmas() is not None) :
            new_array = new_array.select(new_array.sigmas() > 0)
        else :
          raise Sorry("remove_negatives not applicable to %s." % array_name)
      elif (array_params.massage_intensities) :
        if (new_array.is_xray_intensity_array()) :
          if array_params.output_as == "amplitudes" :
            new_array = new_array.enforce_positive_amplitudes()
          else :
            raise Sorry(("You must output %s as amplitudes to use the "+
              "massage_intensities option.") % array_name)
        else :
          raise Sorry("The parameter massage_intensities is only valid for "+
            "X-ray intensity arrays.")
      if (array_params.filter_by_signal_to_noise is not None) :
        if (not is_experimental_data) :
          raise Sorry(("Filtering by signal-to-noise is only supported for "+
            "amplitudes or intensities (failed on %s:%s).") %
            (file_name, array_params.labels))
        elif (array_params.filter_by_signal_to_noise <= 0) :
          raise Sorry(("A value greater than zero is required for the "+
            "cutoff for filtering by signal to noise ratio (failed array: "+
            "%s:%s).") % (file_name, array_params.labels))
        sigmas = new_array.sigmas()
        if (sigmas is None) :
          raise Sorry(("Sigma values must be defined to filter by signal "+
            "to noise ratio (failed on %s:%s).") % (file_name,
              array_params.labels))
        elif (not sigmas.all_ne(0.0)) :
          # XXX should it just remove these too?
          raise Sorry(("The sigma values for the array %s:%s include one or "+
            "more zeros - filtering by signal to noise not supported.") %
            (file_name, array_params.labels))
        data = new_array.data()
        new_array = new_array.select(
          (data / sigmas) > array_params.filter_by_signal_to_noise)
      # leave the default as [0]*6 to make the format clear, but reset to
      # None if unchanged
      if (array_params.add_b_aniso == [0,0,0,0,0,0]) :
        array_params.add_b_aniso = None
      if ((array_params.add_b_iso is not None) or
          (array_params.add_b_aniso is not None)) :
        if (not isinstance(new_array.data(), flex.double) and
                isinstance(new_array.data(), flex.complex_double)) :
          raise Sorry(("Applying a B-factor to the data in %s is not "+
            "permitted.") % array_name)
        if (array_params.add_b_iso is not None) :
          new_array = new_array.apply_debye_waller_factors(
            b_iso=array_params.add_b_iso,
            apply_to_sigmas=True)
        if (array_params.add_b_aniso is not None) :
          new_array = new_array.apply_debye_waller_factors(
            b_cart=array_params.add_b_aniso,
            apply_to_sigmas=True)
      if (array_params.shuffle_values) :
        print >> log, "Shuffling values for %s" % array_name
        combined_array = None
        tmp_array = new_array.deep_copy()
        tmp_array.setup_binner(n_bins=min(100, tmp_array.indices().size()))
        for i_bin in tmp_array.binner().range_used() :
          bin_sel = tmp_array.binner().selection(i_bin)
          bin_array = tmp_array.select(bin_sel)
          perm = flex.random_permutation(bin_array.data().size())
          sigmas = bin_array.sigmas()
          if (sigmas is not None) :
            sigmas = sigmas.select(perm)
          data = bin_array.data().select(perm)
          bin_array = bin_array.customized_copy(data=data, sigmas=sigmas)
          if (combined_array is None) :
            combined_array = bin_array.deep_copy()
          else :
            combined_array = combined_array.concatenate(bin_array)
        if (combined_array.indices().size() != new_array.indices().size()) :
          raise RuntimeError("Array size changed: %d versus %d" %
            (combined_array.indices().size(), new_array.indices().size()))
        new_array = combined_array
      if (array_params.reset_values_to) :
        if (not isinstance(new_array.data(), flex.double) and
            not isinstance(new_array.data(), flex.int)) :
          raise Sorry("Resetting the values for %s is not permitted." %
            array_name)
        print >> log, "Resetting values for %s to %g" % (array_name,
          array_params.reset_values_to)
        data = new_array.data().deep_copy()
        new_value = array_params.reset_values_to
        if isinstance(data, flex.int) :
          new_value = int(new_value)
        data.fill(new_value)
        new_array = new_array.customized_copy(data=data)

      #-----------------------------------------------------------------
      # MISCELLANEOUS
      if new_array.is_xray_intensity_array() :
        if (array_params.output_as in ["amplitudes", "amplitudes_fw"]) :
          if (array_params.output_as == "amplitudes") :
            output_array = new_array.f_sq_as_f()
          else :
            from cctbx import french_wilson
            output_array = french_wilson.french_wilson_scale(
              miller_array=new_array,
              log=log)
          output_array.set_observation_type_xray_amplitude()
          if (array_types[i] is not None) :
            array_types[i] = re.sub("J", "F", array_types[i])
          start_label = None
          if ((output_labels is not None) and
              (output_labels[0].upper().startswith("I"))) :
            raise Sorry(("The output labels for the array %s:%s (%s) are not "+
              "suitable for amplitudes; please change them to something "+
              "with an 'F', or leave this array as intensities") %
              (file_name, array_params.labels, " ".join(output_labels)))
      elif new_array.is_xray_amplitude_array() :
        if array_params.output_as == "intensities" :
          output_array = new_array.f_as_f_sq()
          output_array.set_observation_type_xray_intensity()
          if (array_types[i] is not None) :
            array_types[i] = re.sub("F", "J", array_types[i])
          if ((output_labels is not None) and
              array_params.column_root_label.upper().startswith("F")) :
              raise Sorry("The root label '%s' is not appropriate for "+
                "amplitudes." % array_params.column_root_label)
          elif output_labels[0].upper().startswith("F") :
            raise Sorry(("The output labels for the array %s:%s (%s) are not "+
              "suitable for intensities; please change them to something "+
              "with an 'I', or leave this array as amplitudes.") %
              (file_name, array_params.labels, " ".join(output_labels)))
      if (array_params.force_type != "auto") :
        if (not is_experimental_data) :
          raise Sorry(("You may only override the output observation type for "+
            "amplitudes or intensities - the data in %s:%s are unsupported.") %
            (file_name, array_params.labels))
        if (array_params.force_type == "amplitudes") :
          output_array = new_array.set_observation_type_xray_amplitude()
          if (array_types[i] is not None) :
            array_types[i] = re.sub("J", "F", array_types[i])
        elif (array_params.force_type == "intensities") :
          output_array = new_array.set_observation_type_xray_intensity()
          if (array_types[i] is not None) :
            array_types[i] = re.sub("F", "J", array_types[i])
      if output_array is None :
        output_array = new_array
      if (len(params.mtz_file.exclude_reflection) > 0) :
        for hkl in params.mtz_file.exclude_reflection :
          output_array = output_array.delete_index(hkl)

      # check column_root_label
      if (array_params.column_root_label is not None) :
        root_label = array_params.column_root_label.upper()
        if (output_array.is_xray_intensity_array()) :
          if (not root_label.startswith("I")) :
            raise Sorry(("The specified base column label for the array '%s' "+
              "is inconsistent with the output array type (intensities). "+
              "Please use 'I' (either case) as the first character in the "+
              "label.") % (array_name))
        elif (output_array.is_xray_amplitude_array()) :
          if (not root_label.startswith("F")) :
            raise Sorry(("The specified base column label for the array '%s' "+
              "is inconsistent with the output array type (amplitudes). "+
              "Please use 'F' (either case) as the first character in the "+
              "label.") % (array_name))
        else :
          if (root_label == "I") or (root_label == "F") :
            raise Sorry(("You have specified the base column label '%s' "+
              "for the array '%s', which is neither intensities nor "+
              "amplitudes; the base labels 'I' and 'F' are reserverd "+
              "for these array data types .") %
              (array_params.column_root_label, array_name))

      #-----------------------------------------------------------------
      # OUTPUT
      assert isinstance(output_array, cctbx.miller.array)
      default_label = array_params.column_root_label
      if (output_labels is not None) :
        default_label = 2 * string.uppercase[i]
      column_types = None
      import iotbx.mtz
      default_types = iotbx.mtz.default_column_types(output_array)
      if (array_types[i] is not None) :
        if len(default_types) == len(array_types[i]) :
          print >> log, "Recovering original column types %s" % array_types[i]
          column_types = array_types[i]
      if (output_array.data().size() == 0) :
        raise Sorry("The array %s:%s ended up empty.  Please check the "+
          "resolution cutoffs to make sure they do not exclude all data "+
          "from the input file.  If you think the parameters were correct, "+
          "this is probably a bug; please contact bugs@phenix-online.org "+
          "with a description of the problem.")
      if is_rfree_array(new_array, info) :
        r_free_arrays.append((new_array, info, output_labels,
          array_params.column_root_label, file_name))
      else :
        if DEBUG :
          print >> log, "  Final size:    %d" % output_array.data().size()
          if output_array.sigmas() is not None :
            print >> log, "      sigmas:    %d" % output_array.sigmas().size()
        self.add_array_to_mtz_dataset(
          output_array=output_array,
          default_label=default_label,
          column_types=column_types,
          out=log)
        if (output_labels is not None) :
          for label in output_labels :
            labels.append(label)
            label_files.append(file_name)
        else :
          n_cols = len(default_types)
          if (output_array.anomalous_flag() and
              not output_array.is_xray_reconstructed_amplitude_array()) :
            n_cols *= 2
          for k in range(n_cols) :
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
      if (len(self.final_arrays) > 1) :
        warnings.warn("Multiple Miller arrays are already present in this "+
          "file; the R-free flags will be generated based on the total "+
          "of reflections in all arrays combined.  If you want the fraction "+
          "of test set reflections to be relative to a specific array, you "+
          "should run the editor with that array separately first.",
          UserWarning)
      if (r_free_params.relative_to_complete_set) :
        combined_set = complete_set
      elif (combined_set is not None) :
        self.check_and_warn_about_incomplete_r_free_flags(combined_set)
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
        if (test_flag_value is None) :
          if (r_free_params.old_test_flag_value is not None) :
            test_flag_value = r_free_params.old_test_flag_value
          elif ((r_free_params.warn_if_all_same_value) or
                (r_free_params.extend)) :
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
        if (r_free_params.preserve_input_values) :
          r_free_flags = new_array
        else :
          new_data = (new_array.data()==test_flag_value)
          assert isinstance(new_data, flex.bool)
          r_free_flags = new_array.array(data=new_data)
        r_free_flags = r_free_flags.map_to_asu()
        generate_bijvoet_mates = False
        if not r_free_flags.is_unique_set_under_symmetry() :
          r_free_flags = r_free_flags.merge_equivalents().array()
        if (r_free_flags.anomalous_flag()) :
          r_free_flags = r_free_flags.average_bijvoet_mates()
          if (output_labels is not None) and (len(output_labels) != 1) :
            assert (not combined_set.anomalous_flag())
            # XXX can't do this operation on a miller set - will expand the
            # r-free flags later
            generate_bijvoet_mates = True
        if (r_free_params.adjust_fraction) :
          print >> log, "Resizing test set in %s" % array_name
          r_free_as_bool = get_r_free_as_bool(r_free_flags,
            test_flag_value)
          r_free_flags = r_free_utils.adjust_fraction(
            miller_array=r_free_as_bool,
            fraction=r_free_params.fraction,
            log=log)
        if (r_free_params.extend) :
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
        if (generate_bijvoet_mates) :
          output_array = output_array.generate_bijvoet_mates()
        if (len(params.mtz_file.exclude_reflection) > 0) :
          for hkl in params.mtz_file.exclude_reflection :
            output_array = output_array.delete_index(hkl)
        if (r_free_params.export_for_ccp4) :
          print >> log, "%s: converting to CCP4 convention" % array_name
          output_array = export_r_free_flags(
            miller_array=output_array,
            test_flag_value=True)
        default_label = root_label
        if (default_label is None) :
          default_label = "A" + string.uppercase[i+1]
        self.add_array_to_mtz_dataset(
          output_array=output_array,
          default_label=default_label,
          column_types="I",
          out=log)
        if (output_labels is not None) :
          for label in output_labels :
            labels.append(label)
            label_files.append(file_name)
        else :
          if (output_array.anomalous_flag()) :
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
        r_free_params.force_generate) :
      if (len(self.final_arrays) > 1) :
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
      if (r_free_params.relative_to_complete_set) :
        combined_set = complete_set
      else :
        self.check_and_warn_about_incomplete_r_free_flags(combined_set)
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
        print >> log, "%s: converting to CCP4 convention" % array_name
        output_array = export_r_free_flags(
          miller_array=new_r_free_array,
          test_flag_value=True)
      else:
        output_array = new_r_free_array
      if (len(params.mtz_file.exclude_reflection) > 0) :
        for hkl in params.mtz_file.exclude_reflection :
          output_array = output_array.delete_index(hkl)
      self.add_array_to_mtz_dataset(
        output_array=output_array,
        default_label=r_free_params.new_label,
        column_types="I",
        out=log)
      labels.append(None)
      label_files.append("(new array)")
      self.created_r_free = True

    #-------------------------------------------------------------------
    # RE-LABEL COLUMNS
    mtz_object = self.mtz_dataset.mtz_object()
    self.label_changes = []
    self.mtz_object = mtz_object
    if not len(labels) == mtz_object.n_columns() :
      print >> log, "\n".join([ "LABEL: %s" % label for label in labels ])
      self.show(out=log)
      raise Sorry("The number of output labels does not match the final "+
        "number of labels in the MTZ file.  Details have been printed to the "+
        "console.")
    i = 0
    used = dict([ (label, 0) for label in labels ])
    invalid_chars = re.compile("[^A-Za-z0-9_\-+\(\)]")
    for column in self.mtz_object.columns() :
      if (labels[i] is not None) and (column.label() != labels[i]) :
        label = labels[i]
        original_label = label
        if invalid_chars.search(label) is not None :
          raise Sorry(("Invalid label '%s'.  Output labels may only contain "+
            "alphanumeric characters, underscore, plus and minus signs, or "+
            "parentheses.")
            % label)
        if used[label] > 0 :
          if params.mtz_file.resolve_label_conflicts :
            if label.endswith("(+)") or label.endswith("(-)") :
              label = label[0:-3] + ("_%d" % (used[label]+1)) + label[-3:]
            else :
              label += "_%d" % (used[labels[i]] + 1)
            self.label_changes.append((label_files[i], labels[i], label))
          else :
            raise Sorry(("Duplicate column label '%s'.  Specify "+
              "resolve_label_conflicts=True to automatically generate "+
              "non-redundant labels, or edit the parameters to provide your"+
              "own choice of labels.") % labels[i])
        try :
          column.set_label(label)
        except RuntimeError, e :
          if ("new_label is used already" in str(e)) :
            col_names = [ col.label() for col in mtz_object.columns() ]
            raise RuntimeError(("Duplicate column label '%s': current labels "+
              "are %s; user-specified output labels are %s.") %
              (label, " ".join(col_names), " ".join(labels)))
        else :
          used[original_label] += 1
      i += 1

  def add_array_to_mtz_dataset (self, output_array, default_label,
      column_types, out=sys.stdout) :
    # apply change of basis here
    if (self.params.mtz_file.crystal_symmetry.change_of_basis is not None) :
      from iotbx.reflection_file_converter import apply_change_of_basis
      output_array, cb_op = apply_change_of_basis(
        miller_array=output_array,
        change_of_basis=self.params.mtz_file.crystal_symmetry.change_of_basis,
        eliminate_invalid_indices=\
          self.params.mtz_file.crystal_symmetry.eliminate_invalid_indices,
        out=out)
    if (self.params.mtz_file.crystal_symmetry.eliminate_sys_absent) :
      output_array = output_array.eliminate_sys_absent()
    if self.mtz_dataset is None :
      self.mtz_dataset = output_array.as_mtz_dataset(
        column_root_label=default_label,
        column_types=column_types)
    else :
     self.mtz_dataset.add_miller_array(
        miller_array=output_array,
        column_root_label=default_label,
        column_types=column_types)

  def show (self, out=sys.stdout) :
    if self.mtz_object is not None :
      print >> out, ""
      print >> out, ("=" * 20) + " Summary of output file " + ("=" * 20)
      self.mtz_object.show_summary(out=out, prefix="  ")
      print >> out, ""

  def finish (self) :
    assert self.mtz_object is not None
    if self.params.verbose :
      self.show(out=self.log)
    self.mtz_object.write(file_name=self.params.mtz_file.output_file)
    n_refl = self.mtz_object.n_reflections()
    del self.mtz_object
    self.mtz_object = None
    return n_refl

  def check_and_warn_about_incomplete_r_free_flags (self, combined_set) :
    completeness = combined_set.completeness()
    if (completeness < 0.99) :
      warnings.warn(("The arrays in the input file are incomplete "+
        "(%.1f%% of reflections missing), so the newly generated R-free "+
        "flags will be incomplete as well.  This will not cause a problem "+
        "for refinement, but may result in inconsistent flags later on if "+
        "you collect additional data.  You can disable this warning by "+
        "clicking the box labeled \"Generate flags relative to maximum "+
        "complete set\" in the R-free options dialog window.") %
        (100 * (1-completeness)), UserWarning)

#-----------------------------------------------------------------------
# TODO get rid of these two (need to make sure they aren't imported elsewhere)
def get_r_free_stats (*args, **kwds) :
  from cctbx import r_free_utils
  return r_free_utils.get_r_free_stats(*args, **kwds)

def get_r_free_as_bool (*args, **kwds) :
  from cctbx import r_free_utils
  return r_free_utils.get_r_free_as_bool(*args, **kwds)

def get_best_resolution (miller_arrays, input_symm=None) :
  best_d_min = None
  best_d_max = None
  for array in miller_arrays :
    array_symm = array.crystal_symmetry()
    if ((None in [array_symm.space_group(), array_symm.unit_cell()]) and
        (input_symm is not None)) :
      array = array.customized_copy(crystal_symmetry=input_symm)
    try :
      (d_max, d_min) = array.d_max_min()
      if best_d_max is None or d_max > best_d_max :
        best_d_max = d_max
      if best_d_min is None or d_min < best_d_min :
        best_d_min = d_min
    except Exception, e :
      pass
  return (best_d_max, best_d_min)

def is_rfree_array (miller_array, array_info) :
  from iotbx import reflection_file_utils
  return ((miller_array.is_integer_array() or
           miller_array.is_bool_array()) and
          reflection_file_utils.looks_like_r_free_flags_info(array_info))

def export_r_free_flags (miller_array, test_flag_value) :
  from cctbx import r_free_utils
  new_flags = r_free_utils.export_r_free_flags_for_ccp4(
    flags=miller_array.data(),
    test_flag_value=test_flag_value)
  return miller_array.customized_copy(data=new_flags)

def get_original_array_types (input_file, original_labels) :
  array_types = ""
  mtz_file = input_file.file_object.file_content()
  mtz_columns = mtz_file.column_labels()
  mtz_types = mtz_file.column_types()
  mtz_crossref = dict(zip(mtz_columns, mtz_types))
  for label in original_labels :
    array_types += mtz_crossref[label]
  return array_types

def guess_array_output_labels (miller_array) :
  info = miller_array.info()
  assert info is not None
  labels = info.labels
  output_labels = labels
  if (labels in [["i_obs","sigma"], ["Intensity+-","SigmaI+-"]]) :
    if miller_array.anomalous_flag() :
      output_labels = ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    else :
      output_labels = ["I", "SIGI"]
  elif (miller_array.is_xray_reconstructed_amplitude_array()) :
    output_labels = ["F", "SIGF", "DANO", "SIGDANO", "ISYM"]
  elif ((miller_array.is_xray_amplitude_array() or
         miller_array.is_xray_intensity_array()) and
        miller_array.anomalous_flag()) :
    if (len(labels) == 2) and (miller_array.sigmas() is not None) :
      output_labels = [ "%s(+)" % labels[0],  # data(+)
                        "%s(+)" % labels[1],  # sigma(+)
                        "%s(-)" % labels[0],  # data(-)
                        "%s(-)" % labels[1] ] # sigma(-)
    elif (len(labels) == 1) and (miller_array.sigmas() is None) :
      output_labels = [ "%s(+)" % labels[0], "%s(-)" % labels[0] ]
  return output_labels

# XXX the requirement for defined crystal symmetry in phil input files is
# problematic for automation, where labels and operations are very
# standardized.  so I added this to identify unique symmetry information
# in the collected input files.
def collect_symmetries (file_names) :
  from iotbx import crystal_symmetry_from_any
  file_names = set(file_names)
  file_symmetries = []
  for file_name in file_names :
    symm = crystal_symmetry_from_any.extract_from(file_name)
    if (symm is not None) :
      file_symmetries.append(symm)
  return file_symmetries

def resolve_symmetry (file_symmetries, current_space_group, current_unit_cell):
  space_groups = []
  unit_cells = []
  for symm in file_symmetries :
    if (symm is not None) :
      space_group = symm.space_group_info()
      unit_cell = symm.unit_cell()
      if (space_group is not None) and (unit_cell is not None) :
        if (current_space_group is None) :
          group = space_group.group()
          for other_sg in space_groups :
            other_group = other_sg.group()
            if (other_sg.type().number() != group.type().number()) :
              raise Sorry("Ambiguous space group information in input files - "+
                "please specify symmetry parameters.")
        if (current_unit_cell is None) :
          for other_uc in unit_cells :
            if (not other_uc.is_similar_to(other_uc, 0.001, 0.1)) :
              raise Sorry("Ambiguous unit cell information in input files - "+
                "please specify symmetry parameters.")
        space_groups.append(space_group)
        unit_cells.append(unit_cell)
  consensus_space_group = current_space_group
  consensus_unit_cell = current_unit_cell
  if (len(space_groups) > 0) and (current_space_group is None) :
    consensus_space_group = space_groups[0]
  if (len(unit_cells) > 0) and (current_unit_cell is None) :
    consensus_unit_cell = unit_cells[0]
  return (consensus_space_group, consensus_unit_cell)

def usage (out=sys.stdout, attributes_level=0) :
  print >> out, """
# usage: iotbx.reflection_file_editor [file1.mtz ...] [parameters.eff]
#            --help      (print this message)
#            --details   (show parameter help strings)
# Dumping default parameters:
"""
  master_phil.show(out=out, attributes_level=attributes_level)

def generate_params (file_name, miller_array, include_resolution=False) :
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

def validate_params (params) :
  if (len(params.mtz_file.miller_array) == 0) :
    raise Sorry("No Miller arrays have been selected for the output file.")
  elif len(params.mtz_file.miller_array) > 25 :
    raise Sorry("Only 25 or fewer arrays may be used.")
  if None in [params.mtz_file.crystal_symmetry.space_group,
              params.mtz_file.crystal_symmetry.unit_cell] :
    raise Sorry("Missing or incomplete symmetry information.")
  if (params.mtz_file.r_free_flags.export_for_ccp4 and
      params.mtz_file.r_free_flags.preserve_input_values) :
    raise Sorry("r_free_flags.preserve_input_values and "+
      "r_free_flags.export_for_ccp4 may not be used together.")
  if (params.mtz_file.r_free_flags.adjust_fraction and
      params.mtz_file.r_free_flags.preserve_input_values) :
    raise Sorry("Preserving input values of R-free flags is not supported "+
      "when resizing a test set to the specified fraction.")
  check_if_output_directory_exists(file_name=params.mtz_file.output_file)

#-----------------------------------------------------------------------
def run (args, out=sys.stdout) :
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
    elif os.path.isfile(arg) :
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
      if arg.startswith("--") :
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
      if (not file_name in reflection_file_names) :
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
    for file_name in os.listdir(os.getcwd()) :
      if file_name.startswith("reflections_") and file_name.endswith(".mtz") :
        n += 1
    params.mtz_file.output_file = "reflections_%d.mtz" % (n+1)
  params.mtz_file.output_file = os.path.abspath(params.mtz_file.output_file)
  process = process_arrays(params, input_files=input_files, log=out)
  if params.dry_run :
    print >> out, "# showing final parameters"
    master_phil.format(python_object=params).show(out=out)
    if params.verbose :
      process.show(out=out)
    return process
  process.finish()
  print >> out, "Data written to %s" % params.mtz_file.output_file
  return process

if __name__ == "__main__" :
  run(sys.argv[1:])

#---end
