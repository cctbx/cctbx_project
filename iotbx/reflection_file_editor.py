from __future__ import division

# TODO: regression testing
# TODO: confirm old_test_flag_value if ambiguous

from iotbx import reflection_file_reader, reflection_file_utils, file_reader
from iotbx.reflection_file_utils import get_r_free_flags_scores
import iotbx.phil
import iotbx.mtz
from cctbx import crystal, miller, sgtbx
from scitbx.array_family import flex
from libtbx.phil.command_line import argument_interpreter
from libtbx.math_utils import iceil
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args
import sys, os, string, re, random

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
      .style = bold
    output_space_group = None
      .type = space_group
      .style = bold
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
  resolve_label_conflicts = False
    .type = bool
    .help = Updates label names as necessary to avoid conflicts
    .short_caption = Automatically resolve output label conflicts
    .style = noauto bold
  miller_array
    .multiple = True
    .short_caption = Output Miller array
    .style = fixed auto_align
  {
    file_name = None
      .type = path
      .style = noedit bold
    labels = None
      .type = str
      .short_caption = Array name
      .style = noedit bold
    d_min = None
      .type = float
      .short_caption = High resolution
    d_max = None
      .type = float
      .short_caption = Low resolution
    output_as = *auto amplitudes intensities
      .type = choice
      .short_caption = Output diffraction data as
      .help = If the Miller array is amplitudes or intensities, this flag \
        determines the output data type.
    scale_max = None
      .type = float
      .short_caption = Scale to maximum value
      .help = Scales data such that the maximum is equal to the given value
    scale_factor = None
      .type = float
      .help = Multiplies data with the given factor
    remove_negatives = False
      .type = bool
      .short_caption = Remove negative intensities
    massage_intensities = False
      .type = bool
    output_non_anomalous = False
      .type = bool
      .short_caption = Output non-anomalous data
      .help = If enabled, anomalous arrays will be merged first.  Note that \
        this will cut the number of output labels in half.
    output_label = None
      .type = str
      .multiple = True
      .optional = True
      .short_caption = Output column label
      .help = Most Miller arrays have more than one label, and there must be \
              exactly as many new labels as the number of labels in the \
              old array.  Note however that the output labels do not \
              necessarily correspond to the original array name.  (See caveat \
              in Phenix manual about Scalepack files.)
      .style = fixed
  }
  r_free_flags
    .short_caption = R-free flags generation
    .style = menu_item auto_align box
  {
    generate = True
      .type = bool
      .short_caption = Generate R-free flags if not already present
      .style = bold
    force_generate = False
      .type = bool
      .short_caption = Generate R-free flags even if they are already present
    new_label = FreeR_flag
      .type = str
      .short_caption = Output label for new R-free flags
      .style = bold
    include scope cctbx.miller.generate_r_free_params_str
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
    preserve_input_values = False
      .type = bool
      .short_caption = Preserve original flag values
      .help = If True, CCP4-style R-free flags will be left as random \
        integers instead of being converted to a boolean array.  This option \
        is not compatible with the 'extend' option.
  }
}""", process_includes=True)

class process_arrays (object) :
  def __init__ (self, params, input_files=None, log=sys.stderr,
      accumulation_callback=None, symmetry_callback=None) :
    adopt_init_args(self, locals())
    if len(params.mtz_file.miller_array) == 0 :
      raise Sorry("No Miller arrays have been selected for the output file.")
    elif len(params.mtz_file.miller_array) > 25 :
      raise Sorry("Only 25 or fewer arrays may be used.")
    if None in [params.mtz_file.crystal_symmetry.space_group,
                params.mtz_file.crystal_symmetry.unit_cell] :
      raise Sorry("Missing or incomplete symmetry information.")
    if (params.mtz_file.r_free_flags.extend and
        params.mtz_file.r_free_flags.preserve_input_values) :
      raise Sorry("r_free_flags.preserve_input_values and r_free_flags.extend"+
        " may not be used together.")

    #-------------------------------------------------------------------
    # COLLECT ARRAYS
    miller_arrays = []
    file_names = []
    array_types = []
    for array_params in params.mtz_file.miller_array :
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
    (d_max, d_min) = get_best_resolution(miller_arrays)
    if d_max is None and params.mtz_file.d_max is None :
      raise Sorry("No low-resolution cutoff could be found in the "+
        "parameters or input file(s); you need to explicitly set this value "+
        "for the program to run.")
    if d_min is None :
      if params.mtz_file.d_min is None :
        raise Sorry("No high-resolution cutoff could be found in the "+
          "parameters or input file(s); you need to explicitly set this "+
          "value for the program to run.")
    if params.mtz_file.d_max is not None and params.mtz_file.d_max < d_max :
      d_max = params.mtz_file.d_max
    if params.mtz_file.d_min is not None and params.mtz_file.d_min > d_min :
      d_min = params.mtz_file.d_min
    if params.mtz_file.r_free_flags.random_seed is not None :
      random.seed(params.mtz_file.r_free_flags.random_seed)
      flex.set_random_seed(params.mtz_file.r_free_flags.random_seed)

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
      space_group=params.mtz_file.crystal_symmetry.space_group.group())
    derived_sg = input_symm.space_group().build_derived_point_group()
    output_symm = crystal.symmetry(unit_cell=output_uc,
                                   space_group=output_sg)

    #-------------------------------------------------------------------
    # MAIN LOOP
    i = 0
    special_labels = ["i_obs,sigma", "Intensity+-,SigmaI+-"]
    r_free_arrays = []
    for (array_params, file_name, miller_array) in \
        zip(params.mtz_file.miller_array, file_names, miller_arrays) :
      array_name = "%s:%s" % (file_name, array_params.labels)
      if params.verbose :
        print "Processing %s" % array_name
      if array_params.d_max is not None and array_params.d_max <= 0 :
        array_params.d_max = None
      if array_params.d_min is not None and array_params.d_min <= 0 :
        array_params.d_min = None
      output_array = None # this will eventually be the final processed array
      output_labels = array_params.output_label
      info = miller_array.info()
      if not None in [array_params.scale_factor, array_params.scale_max] :
        raise Sorry("The parameters scale_factor and scale_max are " +
          "mutually exclusive.")
      if not False in [array_params.remove_negatives,
                       array_params.massage_intensities] :
        raise Sorry("The parameters remove_negatives and massage_intensities "+
          "are mutually exclusive.")
      if DEBUG :
        print "  Starting size:  %d" % miller_array.data().size()
        if miller_array.sigmas() is not None :
          print "         sigmas:  %d" % miller_array.sigmas().size()

      #-----------------------------------------------------------------
      # OUTPUT LABELS SANITY CHECK
      if miller_array.anomalous_flag() :
        labels_base = re.sub(",merged$", "", array_params.labels)
        input_labels = labels_base.split(",")
        if array_params.output_non_anomalous :
          if labels_base in special_labels and not len(output_labels) == 2 :
            raise Sorry(("There are too many output labels for the array "+
              "%s, which is being converted to non-anomalous data. "+
              "Labels such as I,SIGI are appropriate (or F,SIGF if you are "+
              "converting the array to amplitudes.") % array_name)
          elif len(output_labels) == len(input_labels) :
            raise Sorry(("There are too many output labels for the array "
              "%s, which is being converted to non-anomalous data. "+
              "The total number of columns will be halved in the output "+
              "array, and the labels should not have trailing (+) or (-).") %
              array_name)
        elif labels_base in special_labels and not len(output_labels) == 4 :
          raise Sorry(("There are not enough output labels for the array "+
            "%s. For Scalepack or d*TREK files containing anomalous "+
            "data, you must specify exactly four column labels (e.g. "+
            "I(+),SIGI(+),I(-),SIGI(-), or F(+),SIGF(+),F(-),SIGF(-) "+
            "if you are converting the array to amplitudes).") %
            array_name)

      #-----------------------------------------------------------------
      # APPLY SYMMETRY
      array_sg = miller_array.space_group()
      array_uc = miller_array.unit_cell()
      if array_sg is not None :
        if array_sg.build_derived_point_group() != derived_sg :
          raise Sorry(("The point group for the Miller array %s (%s) does "+
            "not match the point group of the overall space group (%s).") %
            (array_name, str(array_sg), str(input_symm.space_group())))
      if array_uc is not None :
        if not array_uc.is_similar_to(input_symm.unit_cell()) :
          raise Sorry(("The unit cell for the Miller array %s (%s) is "+
            "significantly different than the output unit cell (%s).") %
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
          raise Sorry(("The data in %s cannot be merged because it is in "+
            "integer format.  If you wish to change symmetry (or the input "+
            "data is unmerged), you must omit this array.  (Note also that "+
            "merging will fail for R-free flags if the flags for symmetry-"+
            "related reflections are not identical.)") % array_name)
        new_array = new_array.merge_equivalents().array()

      if DEBUG :
        print "  Adjusted size:  %d" % new_array.data().size()
        if miller_array.sigmas() is not None :
          print "         sigmas:  %d" % new_array.sigmas().size()

      #-----------------------------------------------------------------
      # CHANGE OF BASIS
      # XXX: copied from reflection_file_converter nearly verbatim
      if params.mtz_file.crystal_symmetry.change_of_basis is not None :
        if change_symmetry :
          raise Sorry("You may not change symmetry when change_of_basis is "+
            "defined.")
        c_o_b = params.mtz_file.crystal_symmetry.change_of_basis
        if c_o_b == "to_reference_setting" :
          cb_op = new_array.change_of_basis_op_to_reference_setting()
        elif c_o_b == "to_primitive_setting" :
          cb_op = new_array.change_of_basis_op_to_primitive_setting()
        elif c_o_b == "to_niggli_cell":
          cb_op = new_array.change_of_basis_op_to_niggli_cell()
        elif c_o_b == "to_inverse_hand" :
          cb_op = new_array.change_of_basis_op_to_inverse_hand()
        else:
          cb_op = sgtbx.change_of_basis_op(c_o_b)
        if (cb_op.c_inv().t().is_zero()):
          print >> log, ("  Change of basis operator in both h,k,l and "+
                         "x,y,z notation:")
          print >> log, "   ", cb_op.as_hkl()
        else:
          print >> log, "  Change of basis operator in x,y,z notation:"
        print >> log, "    %s [Inverse: %s]" % (cb_op.as_xyz(),
          cb_op.inverse().as_xyz())
        if (d < 0 and co.change_of_basis != "to_inverse_hand"):
          print >> out, ("WARNING: This change of basis operator changes the "+
                        "hand!")
        if params.mtz_file.crystal_symmetry.eliminate_invalid_indices :
          sel = cb_op.apply_results_in_non_integral_indices(
            miller_indices=new_array.indices())
          toss = flex.bool(new_array.indices().size(),sel)
          keep = ~toss
          keep_array = new_array.select(keep)
          toss_array = new_array.select(toss)
          print >> out, "  Mean value for kept reflections:", \
            flex.mean(keep_array.data())
          print >> out, "  Mean value for invalid reflections:", \
            flex.mean(toss_array.data())
          new_array = new_array
        new_array = new_array.change_basis(cb_op=cb_op)
        print >> out, "  Crystal symmetry after change of basis:"
        crystal.symmetry.show_summary(new_array, out=out, prefix="    ")

      #-----------------------------------------------------------------
      # OTHER FILTERING
      if DEBUG :
        print "  Resolution before array-specific filter: %.2f - %.2f" % (
          new_array.d_max_min())
      # first the array-specific cutoff
      new_array = new_array.resolution_filter(
        d_min=array_params.d_min,
        d_max=array_params.d_max)
      if DEBUG :
        print "              after array-specific filter: %.2f - %.2f" % (
          new_array.d_max_min())
      # now apply the global cutoff
      new_array = new_array.resolution_filter(
          d_min=params.mtz_file.d_min,
          d_max=params.mtz_file.d_max)
      if DEBUG :
        print "  Truncated size: %d" % new_array.data().size()
        if new_array.sigmas() is not None :
          print "          sigmas: %d" % new_array.sigmas().size()
      if new_array.anomalous_flag() and array_params.output_non_anomalous :
        print >> log, ("Converting array %s from anomalous to non-anomalous." %
                       array_name)
        if not new_array.is_xray_intensity_array() :
          new_array = new_array.average_bijvoet_mates()
        else :
          new_array = new_array.f_sq_as_f()
          new_array = new_array.average_bijvoet_mateS()
          new_array = new_array.f_as_f_sq()
          new_array.set_observation_type_xray_intensity()
      if array_params.scale_max is not None :
        print >> log, ("Scaling %s such that the maximum value is: %.6g" %
                       (array_name, array_params.scale_max))
        new_array = new_array.apply_scaling(target_max=array_params.scale_max)
      elif array_params.scale_factor is not None :
        print >> log, ("Multiplying data in %s with the factor: %.6g" %
                       (array_name, array_params.scale_factor))
        new_array = new_array.apply_scaling(factor=array_params.scale_factor)
      if array_params.remove_negatives :
        if new_array.is_real_array() :
          print >> log, "Removing negatives from %s" % array_name
          new_array = new_array.select(new_array.data() > 0)
          if new_array.sigmas() is not None :
            new_array = new_array.select(new_array.sigmas() > 0)
        else :
          raise Sorry("remove_negatives not applicable to %s." % array_name)
      elif array_params.massage_intensities :
        if new_array.is_xray_intensity_array() :
          if array_params.output_as == "amplitudes" :
            new_array = new_array.enforce_positive_amplitudes()
          else :
            raise Sorry(("You must output %s as amplitudes to use the "+
              "massage_intensities option.") % array_name)
        else :
          raise Sorry("The parameter massage_intensities is only valid for "+
            "X-ray intensity arrays.")

      #-----------------------------------------------------------------
      # MISCELLANEOUS
      if new_array.is_xray_intensity_array() :
        if array_params.output_as == "amplitudes" :
          output_array = new_array.f_sq_as_f()
          if output_labels[0].upper().startswith("I") :
            raise Sorry(("The output labels for the array %s:%s (%s) are not "+
              "suitable for amplitudes; please change them to something "+
              "with an 'F', or leave this array as intensities") %
              (file_name, array_params.labels, " ".join(output_labels)))
      elif new_array.is_xray_amplitude_array() :
        if array_params.output_as == "intensities" :
          output_array = new_array.f_as_f_sq()
          if output_labels[0].upper().startswith("F") :
            raise Sorry(("The output labels for the array %s:%s (%s) are not "+
              "suitable for intensities; please change them to something "+
              "with an 'I', or leave this array as amplitudes.") %
              (file_name, array_params.labels, " ".join(output_labels)))
      if output_array is None :
        output_array = new_array

      #-----------------------------------------------------------------
      # OUTPUT
      fake_label = 2 * string.uppercase[i]
      column_types = None
      if array_types[i] is not None :
        default_types = iotbx.mtz.default_column_types(output_array)
        if len(default_types) == len(array_types[i]) :
          print >> log, "Recovering original column types %s" % array_types[i]
          column_types = array_types[i]
      if output_array.data().size() == 0 :
        raise Sorry("The array %s:%s ended up empty.  Please check the "+
          "resolution cutoffs to make sure they do not exclude all data "+
          "from the input file.  If you think the parameters were correct, "+
          "this is probably a bug; please contact bugs@phenix-online.org "+
          "with a description of the problem.")
      if is_rfree_array(new_array, info) :
        r_free_arrays.append((new_array, info, output_labels, file_name))
      else :
        if DEBUG :
          print "  Final size:    %d" % output_array.data().size()
          if output_array.sigmas() is not None :
            print "      sigmas:    %d" % output_array.sigmas().size()
        self.add_array_to_mtz_dataset(
          output_array=output_array,
          fake_label=fake_label,
          column_types=column_types)
        for label in output_labels :
          labels.append(label)
          label_files.append(file_name)
        self.final_arrays.append(output_array)
      i += 1

    #-------------------------------------------------------------------
    # EXISTING R-FREE ARRAYS
    if len(r_free_arrays) > 0 :
      have_r_free_array = True
      if len(self.final_arrays) > 0 :
        complete_set = make_joined_set(self.final_arrays).complete_set()
      else :
        complete_set = None
      i = 0
      for (new_array, info, output_labels, file_name) in r_free_arrays :
        flag_scores = get_r_free_flags_scores(miller_arrays=[new_array],
           test_flag_value=params.mtz_file.r_free_flags.old_test_flag_value)
        test_flag_value = flag_scores.test_flag_values[0]
        if params.mtz_file.r_free_flags.preserve_input_values :
          r_free_flags = new_array
        else :
          r_free_flags = new_array.array(data=new_array.data()==test_flag_value)
        fraction_free = (r_free_flags.data().count(True) /
                         r_free_flags.data().size())
        print >>log, "%s: fraction_free=%.3f" % (info.labels[0], fraction_free)
        if complete_set is not None :
          missing_set = complete_set.lone_set(r_free_flags.map_to_asu())
        else :
          missing_set = r_free_flags.complete_set(d_min=d_min,
            d_max=d_max).lone_set(r_free_flags.map_to_asu())
        n_missing = missing_set.indices().size()
        print >>log, "%s: missing %d reflections" % (info.labels[0], n_missing)
        if n_missing != 0 and params.mtz_file.r_free_flags.extend :
          if n_missing <= 20 :
            # FIXME: MASSIVE CHEAT necessary for tiny sets
            missing_flags = missing_set.array(data=flex.bool(n_missing,False))
          else :
            if accumulation_callback is not None :
              if not accumulation_callback(miller_array=new_array,
                                           test_flag_value=test_flag_value,
                                           n_missing=n_missing,
                                           column_label=info.labels[0]) :
                continue
            missing_flags = missing_set.generate_r_free_flags(
              fraction=fraction_free,
              max_free=None,
              use_lattice_symmetry=True)
          output_array = r_free_flags.concatenate(other=missing_flags)
          if not output_array.is_unique_set_under_symmetry() :
            output_array = output_array.merge_equivalents().array()
        else :
          output_array = r_free_flags
        if params.mtz_file.r_free_flags.export_for_ccp4 :
          print >> log, "%s: converting to CCP4 convention" % array_name
          output_array = export_r_free_flags(miller_array=output_array,
            test_flag_value=True)
        fake_label = "A" + string.uppercase[i+1]
        self.add_array_to_mtz_dataset(
          output_array=output_array,
          fake_label=fake_label,
          column_types="I")
        for label in output_labels :
          labels.append(label)
          label_files.append(file_name)
        self.final_arrays.append(output_array)
        i += 1

    #-------------------------------------------------------------------
    # NEW R-FREE ARRAY
    if ((params.mtz_file.r_free_flags.generate and not have_r_free_array) or
        params.mtz_file.r_free_flags.force_generate) :
      complete_set = make_joined_set(self.final_arrays).complete_set()
      r_free_params = params.mtz_file.r_free_flags
      new_r_free_array = complete_set.generate_r_free_flags(
        fraction=r_free_params.fraction,
        max_free=r_free_params.max_free,
        lattice_symmetry_max_delta=r_free_params.lattice_symmetry_max_delta,
        use_lattice_symmetry=r_free_params.use_lattice_symmetry,
        use_dataman_shells=r_free_params.use_dataman_shells,
        n_shells=r_free_params.n_shells)
      if r_free_params.new_label is None or r_free_params.new_label == "" :
        r_free_params.new_label = "FreeR_flag"
      if params.mtz_file.r_free_flags.export_for_ccp4 :
        print >> log, "%s: converting to CCP4 convention" % array_name
        output_array = export_r_free_flags(miller_array=new_r_free_array,
          test_flag_value=True)
      self.mtz_dataset.add_miller_array(
        miller_array=new_r_free_array,
        column_root_label=r_free_params.new_label)
      labels.append(r_free_params.new_label)
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
      if column.label() != labels[i] :
        label = labels[i]
        if invalid_chars.search(label) is not None :
          raise Sorry(("Invalid label '%s'.  Output labels may only contain "+
            "alphanumeric characters, underscore, plus and minus signs, or "+
            "parentheses.")
            % label)
        if used[labels[i]] > 0 :
          if params.mtz_file.resolve_label_conflicts :
            if label.endswith("(+)") or label.endswith("(-)") :
              label = label[0:-3] + ("_%d" % (used[labels[i]]+1)) + label[-3:]
            else :
              label += "_%d" % (used[labels[i]] + 1)
            self.label_changes.append((label_files[i], labels[i], label))
          else :
            raise Sorry(("Duplicate column label '%s'.  Specify "+
              "resolve_label_conflicts=True to automatically generate "+
              "non-redundant labels, or edit the parameters to provide your"+
              "own choice of labels.") % labels[i])
        column.set_label(label)
        used[labels[i]] += 1
      i += 1

  def add_array_to_mtz_dataset (self, output_array, fake_label, column_types) :
    if self.mtz_dataset is None :
      self.mtz_dataset = output_array.as_mtz_dataset(
        column_root_label=fake_label,
        column_types=column_types)
    else :
     self.mtz_dataset.add_miller_array(
        miller_array=output_array,
        column_root_label=fake_label,
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
    del self.mtz_object
    self.mtz_object = None

#-----------------------------------------------------------------------
def get_r_free_stats (miller_array, test_flag_value) :
  array = get_r_free_as_bool(miller_array, test_flag_value)
  n_free = array.data().count(True)
  accu =  array.sort(by_value="resolution").r_free_flags_accumulation()
  lr = flex.linear_regression(accu.reflection_counts.as_double(),
                              accu.free_fractions)
  assert lr.is_well_defined()
  slope = lr.slope()
  y_ideal = accu.reflection_counts.as_double() * slope
  sse = 0
  n_bins = 0
  n_ref_last = 0
  sse = flex.sum(flex.pow(y_ideal - accu.free_fractions, 2))
  for x in accu.reflection_counts :
    if x > (n_ref_last + 1) :
      n_bins += 1
    n_ref_last = x
  return (n_bins, n_free, sse, accu)

def get_r_free_as_bool (miller_array, test_flag_value=0) :
  if miller_array.is_bool_array() :
    return miller_array
  else :
    assert miller_array.is_integer_array()
    return miller_array.customized_copy(
      data=miller_array.data() == test_flag_value,
      sigmas=None)

def get_best_resolution (miller_arrays) :
  best_d_min = None
  best_d_max = None
  for array in miller_arrays :
    try :
      (d_max, d_min) = array.d_max_min()
      if best_d_max is None or d_max > best_d_max :
        best_d_max = d_max
      if best_d_min is None or d_min < best_d_min :
        best_d_min = d_min
    except Exception, e :
      pass
  return (best_d_max, best_d_min)

def make_joined_set (miller_arrays) :
  master_set = miller.set(
    crystal_symmetry=miller_arrays[0].crystal_symmetry(),
    indices=miller_arrays[0].indices(),
    anomalous_flag=False)
  master_indices = miller_arrays[0].indices().deep_copy()
  for array in miller_arrays[1:] :
    current_indices = array.indices()
    missing_isel = miller.match_indices(master_indices,
      current_indices).singles(1)
    missing_indices = current_indices.select(missing_isel)
    master_indices.extend(missing_indices)
  master_set = miller.set(
    crystal_symmetry=miller_arrays[0].crystal_symmetry(),
    indices=master_indices,
    anomalous_flag=False)
  return master_set.unique_under_symmetry()

def is_rfree_array (miller_array, array_info) :
  return ((miller_array.is_integer_array() or
           miller_array.is_bool_array()) and
          reflection_file_utils.looks_like_r_free_flags_info(array_info))

def export_r_free_flags (miller_array, test_flag_value) :
  if miller_array.is_bool_array() : # XXX: in practice, this is always true
    test_flag_value = True
  else :
    assert miller_array.is_integer_array()
  data = miller_array.data()
  unique_values = set(data)
  if len(unique_values) > 2 : # XXX: is this safe?
    return miller_array
  new_flags = flex.int(data.size())
  for i in range(data.size()) :
    if data[i] == test_flag_value :
      new_flags[i] = 0
    else :
      new_flags[i] = iceil(random.random() * 19)
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
  if labels in [["i_obs","sigma"], ["Intensity+-","SigmaI+-"]] :
    if miller_array.anomalous_flag() :
      output_labels = ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    else :
      output_labels = ["I", "SIGI"]
  return output_labels

def usage (out=sys.stdout, attributes_level=0) :
  print >> out, """
# usage: iotbx.reflection_file_editor [file1.mtz ...] [parameters.eff]
#            --help      (print this message)
#            --details   (show parameter help strings)
# Dumping default parameters:
"""
  master_phil.show(out=out, attributes_level=attributes_level)

def generate_params (file_name, miller_array) :
  param_str = """mtz_file.miller_array {
  file_name = %s
  labels = %s
""" % (file_name, miller_array.info().label_string())
  output_labels = guess_array_output_labels(miller_array)
  for label in output_labels :
    param_str += "  output_label = %s\n" % label
  try :
    (d_max, d_min) = miller_array.d_max_min()
    param_str += """  d_max = %.5f\n  d_min = %.5f\n""" % (d_max, d_min)
  except Exception :
    pass
  param_str += "}"
  return param_str

#-----------------------------------------------------------------------
def run (args, out=sys.stdout) :
  crystal_symmetry_from_pdb = None
  crystal_symmetries_from_hkl = []
  reflection_files = []
  user_phil = []
  all_arrays = []
  if len(args) == 0 :
    usage()
  interpreter = argument_interpreter(master_phil=master_phil,
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
        miller_arrays = input_file.file_object.as_miller_arrays()
        for array in miller_arrays :
          symm = array.crystal_symmetry()
          if symm is not None :
            crystal_symmetries_from_hkl.append(symm)
            break
        reflection_files.append(input_file)
      elif input_file.file_type == "pdb" :
        symm = input_file.file_object.crystal_symmetry()
        if symm is not None :
          crystal_symmetry_from_pdb = symm
    else :
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
    for i, symm in enumerate(crystal_symmetries_from_hkl) :
      params.mtz_file.crystal_symmetry.space_group = symm.space_group_info()
      params.mtz_file.crystal_symmetry.unit_cell = symm.unit_cell()
      break
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
