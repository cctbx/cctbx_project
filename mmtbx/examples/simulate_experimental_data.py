import math
from cctbx import crystal
from cctbx import sgtbx
from cctbx import xray
from cctbx.array_family import flex
from libtbx.utils import Sorry, multi_out
import iotbx.phil
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import pdb
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx import f_model
from cStringIO import StringIO
import sys, os


class error_swap(object):
  def __init__(self,
               miller_obs,
               miller_calc,
               miller_mock,
               n_reso_bins=25,
               n_e_bins = 20,
               thres=3.0):
    self.miller_obs = miller_obs
    self.miller_calc = miller_calc
    self.miller_mock = miller_mock

    # take a common set to avoid possible problems
    self.miller_calc = self.miller_calc.common_set( self.miller_obs )
    self.miller_mock = self.miller_mock.common_set( self.miller_obs )


    # we need to normalise the data, both fobs and fcalc
    norma_obs_obj = absolute_scaling.kernel_normalisation( self.miller_obs,auto_kernel=True )
    norma_calc_obj = absolute_scaling.kernel_normalisation( self.miller_calc,auto_kernel=True )
    norma_mock_obj = absolute_scaling.kernel_normalisation( self.miller_mock,auto_kernel=True )
    self.norma_obs  = norma_obs_obj.normalised_miller_dev_eps.f_sq_as_f()           # normalized data (dived by eps)
    self.norma_calc = norma_calc_obj.normalised_miller_dev_eps.f_sq_as_f()          # as above, for calculated data
    self.norma_mock = norma_mock_obj.normalised_miller_dev_eps.f_sq_as_f()          # as above, for mock data
    self.norma_obs_const =  norma_obs_obj.normalizer_for_miller_array   # the divisor (no eps)
    self.norma_calc_const = norma_calc_obj.normalizer_for_miller_array  # as above
    self.norma_mock_const = norma_mock_obj.normalizer_for_miller_array  # as above

    self.thres = thres

    self.n_reso_bins = n_reso_bins
    self.n_e_bins = n_e_bins
    # first set up a binner please
    self.miller_obs.setup_binner(n_bins = self.n_reso_bins )
    self.miller_calc.use_binner_of( self.miller_obs )
    self.miller_mock.use_binner_of( self.miller_obs )
    self.norma_obs.use_binner_of( self.miller_obs )
    self.norma_calc.use_binner_of( self.miller_calc )
    self.norma_mock.use_binner_of( self.miller_mock )

    self.new_norma_obs = self.norma_obs.deep_copy().set_observation_type( self.norma_obs )

    self.new_obs = None
    self.swap_it()
    #we have to denormalize the data now
    self.new_obs = self.norma_obs.customized_copy(
      data   = self.new_norma_obs.data()*self.new_norma_obs.epsilons().data().as_double()*flex.sqrt(self.norma_calc_const),
      sigmas = self.new_norma_obs.sigmas()*self.new_norma_obs.epsilons().data().as_double()*flex.sqrt(self.norma_calc_const)
    ).set_observation_type( self.miller_obs )
    # all done

  def swap_it(self):
    # for each bin
    for ibin in self.norma_obs.binner().range_used():
      # select all indices in this bin please
      selection = self.norma_obs.binner().bin_indices() == ibin
      tmp_norm_obs = self.norma_obs.select( selection )
      tmp_norm_calc = self.norma_calc.select( selection )
      tmp_norm_mock = self.norma_mock.select( selection )
      # we now have a set of e values, send both arrays of to another routine
      new_e,new_se = self.do_something_clever( tmp_norm_obs.data(),
                                               tmp_norm_obs.sigmas(),
                                               tmp_norm_calc.data(),
                                               tmp_norm_mock.data() )
      self.new_norma_obs = self.norma_obs.customized_copy(
        data = self.norma_obs.data().set_selected(selection,new_e),
        sigmas = self.norma_obs.sigmas().set_selected(selection,new_se)
      )

  def do_something_clever(self,obs,sobs,calc,mock):
    # first get the sort order
    # sort on the calculated data please
    sort_order = flex.sort_permutation( calc )
    inverse_sort_order = sort_order.inverse_permutation()

    sorted_obs  = obs.select(sort_order)
    sorted_sobs = sobs.select(sort_order)
    sorted_calc = calc.select(sort_order)
    sorted_mock = mock.select(sort_order)

    log_calc = flex.log(sorted_mock)
    deltas   = flex.log(sorted_obs) - flex.log(sorted_calc)

    old_deltas = deltas.deep_copy()

    # make bins on the basis of the order
    bin_size = float(sorted_obs.size())/self.n_e_bins
    bin_size = int(bin_size) + 1
    ebin = flex.int()
    count=0
    for ii in xrange( sorted_obs.size() ):
      if ii%bin_size==0:
        count+=1
      ebin.append( count-1 )

    # the bins have been setup, now we can reorder stuff
    for ibin in xrange(self.n_e_bins):
      this_bin_selection = flex.bool( ebin == ibin )
      tmp_n = (this_bin_selection).count(True)
      permute = flex.sort_permutation( flex.random_double( tmp_n ) )

      #select and swap
      selected_deltas = deltas.select( this_bin_selection )
      selected_deltas = selected_deltas.select( permute )
      selected_sobs   = sorted_sobs.select( this_bin_selection )
      selected_sobs   = selected_sobs.select( permute )


      # we have to make a sanity check so that the selected deltas are not very weerd
      # a safeguard to prevent the introductoin of outliers
      mean_delta = flex.mean( selected_deltas )
      std_delta  = math.sqrt( flex.mean( selected_deltas*selected_deltas ) - mean_delta*mean_delta )
      outliers = flex.bool( flex.abs(selected_deltas-mean_delta)>self.thres*std_delta )
      #print list( flex.abs(selected_deltas-mean_delta)/std_delta )
      #print list( outliers )

      if (outliers).count(True) > 0 :
        non_out_delta   = selected_deltas.select( ~outliers )
        tmp_permut      = flex.sort_permutation( flex.random_double( (~outliers).count(True)  ) )
        tmp_delta       = non_out_delta.select( tmp_permut )
        tmp_delta       = tmp_delta[0:(outliers).count(True)]
        selected_deltas = selected_deltas.set_selected( outliers.iselection(), tmp_delta )


      #set the deltas back please
      deltas = deltas.set_selected(this_bin_selection, selected_deltas)
      sorted_sobs = sorted_sobs.set_selected(this_bin_selection, selected_sobs)

    #the deltas have been swapped, apply things back please
    log_calc = log_calc + deltas
    log_calc = flex.exp(log_calc)

    #now we have to get things back in proper order again thank you
    new_fobs = log_calc.select(inverse_sort_order)
    new_sobs = sorted_sobs.select(inverse_sort_order)
    return new_fobs, new_sobs






def select_crystal_symmetry(
      from_command_line,
      from_parameter_file,
      from_coordinate_files,
      from_reflection_files):
  result = crystal.symmetry(
    unit_cell=None,
    space_group_info=None)
  if (from_command_line is not None):
    result = result.join_symmetry(
      other_symmetry=from_command_line, force=False)
  if (from_parameter_file is not None):
    result = result.join_symmetry(
      other_symmetry=from_parameter_file, force=False)
  if (result.unit_cell() is None):
    for crystal_symmetry in from_reflection_files:
      unit_cell = crystal_symmetry.unit_cell()
      if (unit_cell is not None):
        result = crystal.symmetry(
          unit_cell=unit_cell,
          space_group_info=result.space_group_info(),
          assert_is_compatible_unit_cell=False)
        break
  for crystal_symmetry in from_coordinate_files:
    result = result.join_symmetry(other_symmetry=crystal_symmetry, force=False)
  if (result.space_group_info() is None):
    for crystal_symmetry in from_reflection_files:
      space_group_info = crystal_symmetry.space_group_info()
      if (space_group_info is not None):
        result = crystal.symmetry(
          unit_cell=result.unit_cell(),
          space_group_info=space_group_info,
          assert_is_compatible_unit_cell=False)
        break
  return result



master_params = iotbx.phil.parse("""\
simul_utils{
  input{
    unit_cell=None
    .type=unit_cell
    space_group=None
    .type=space_group
    xray_data{
      file_name=None
      .type=path
      labels=None
      .type=str
    }
    model{
      file_name=None
      .type=path
    }
    mock_model{
      file_name=None
      .type=path
    }
  }
  output{
    logfile=simul.log
    .type=str
    hklout=simul.mtz
    .type=str
  }
}
""")

def print_help():
  print """
mmtbx.simulate_data:

Allows one to quickly simulate data 'experimental' data with similar
Fobs-Fcalc distribution as in the given model/data pair.

The keywords are sumarized below and should be self explanatory:

simul_utils{
  input{
    unit_cell=None
    space_group=None
    xray_data{
      file_name=None
      labels=None
    }
    model{
      file_name=None
    }
    mock_model{
      file_name=None
    }
  }
  output{
    logfile=simul.log
    hklout=simul.mtz
  }
}

The main purpose of this file is to generate data with errors that look real.
This is what is done:
1. The pdb file 'model.file_name' is scaled to the observed data 'xray_data.file_name'
   (Fcalc and Fobs are now available and on the same scale)
2. Structure factors are computed for 'mock_model.file_name' with same bulk solvent
   parameters as 'model.file_name' (call this Fmock)
3a. For each resolution bin, generate about 20 E value bins.
3b. In each E value bin, compute ratio=Fobs/Fcalc (or log Fobs - log Fcalc so you will)
3c. make a random permuation of these ratios log differences (call this array random_ratio)
3d. Fmockobs = F_mock*random_ratio
4.  Write out Fmockobs

The data generated in this manner has similar F/sigF (the sigmas are permuted along with the ratios)
and R value properties.

If no mock model is supplied, the model will be the mock model.
The mock model is supposed to be in the same unit cell/spacegroup (this is enforced).

  """


def simul_utils(args):
  if len(args)==0:
    print_help()
  elif ( "--help" in args ):
    print_help()
  elif ( "--h" in args ):
    print_help()
  elif ("-h" in args ):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = master_params.command_line_argument_interpreter(
      home_scope="map_coefs")

    print >> log, "#phil __OFF__"
    print >> log, "======================"
    print >> log, "          SIMUL       "
    print >> log, "A data simulation tool"
    print >> log, "======================"
    print >> log


    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    """
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    """
    # now get the unit cell from the pdb file

    hkl_xs = crystal_symmetry_from_any.extract_from(
      file_name=params.simul_utils.input.xray_data.file_name)
    pdb_xs = crystal_symmetry_from_any.extract_from(
      file_name=params.simul_utils.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.simul_utils.input.unit_cell,
      space_group_info=params.simul_utils.input.space_group  )


    combined_xs = select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[hkl_xs])

    # inject the unit cell and symmetry in the phil scope please
    params.simul_utils.input.unit_cell = combined_xs.unit_cell()
    params.simul_utils.input.space_group = \
      sgtbx.space_group_info( group = combined_xs.space_group() )

    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    print >> log, "#phil __END__"

    if params.simul_utils.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.simul_utils.input.space_group is None:
      raise Sorry("space group not specified")
    if params.simul_utils.input.xray_data.file_name is None:
      raise Sorry("Xray data not specified")
    if params.simul_utils.input.model.file_name is None:
      raise Sorry("pdb file with  model not specified")

    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #
    phil_xs = crystal.symmetry(
      unit_cell=params.simul_utils.input.unit_cell,
      space_group_info=params.simul_utils.input.space_group  )

    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = phil_xs,
      force_symmetry = True,
      reflection_files=[])

    miller_array = None

    miller_array = xray_data_server.get_xray_data(
      file_name = params.simul_utils.input.xray_data.file_name,
      labels = params.simul_utils.input.xray_data.labels,
      ignore_all_zeros = True,
      parameter_scope = 'simul_utils.input.xray_data',
      parameter_name = 'labels'
      )

    info = miller_array.info()

    miller_array = miller_array.map_to_asu()

    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))

    miller_array = miller_array.select(
      miller_array.data() > 0 )
    if  miller_array.sigmas() is not None:
      miller_array = miller_array.select(
        miller_array.sigmas() > 0 )

    if (miller_array.is_xray_intensity_array()):
      miller_array = miller_array.f_sq_as_f()
    elif (miller_array.is_complex_array()):
      miller_array = abs(miller_array)

    miller_array.set_info(info)
    print >> log
    print >> log, "Summary info of observed data"
    print >> log, "============================="
    miller_array.show_summary(f=log)
    print >> log


    free_flags = miller_array.generate_r_free_flags()



    #--------------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    model = pdb.input(file_name=params.simul_utils.input.model.file_name).xray_structure_simple(
                      crystal_symmetry=phil_xs,
                      )
    print >> log, "Atomic model summary"
    print >> log, "===================="
    model.show_summary()
    print >> log

    #-------------------------------------------------------------------
    # Step 3: make an F_model object to get model phases and amplitudes
    #
    print >> log, "Performing bulk solvent scaling"
    print >> log, "==============================="
    print >> log
    print >> log

    f_model_object = f_model.manager(
        f_obs = miller_array,
        r_free_flags = free_flags,
        xray_structure = model )
    f_model_object.update_solvent_and_scale(out=log)
    fmodel = abs( f_model_object.f_model() ).set_observation_type( miller_array )

    mockfmodel = None
    if params.simul_utils.input.mock_model.file_name is not None:
      print >> log, "Reading in mock model"
      print >> log, "====================="
      print >> log
      print >> log
      mock_model = pdb.input(file_name=params.simul_utils.input.mock_model.file_name).xray_structure_simple(
                             crystal_symmetry=phil_xs)
      mock_f_model = f_model.manager(
        f_obs = miller_array,
        r_free_flags = free_flags,
        xray_structure = mock_model )

      mock_f_model.update(
        k_sol  = f_model_object.k_sol() ,
        b_sol  = f_model_object.b_sol() ,
        b_cart = f_model_object.b_cart()
      )
      mockfmodel = abs( mock_f_model.f_model() ).set_observation_type( miller_array )
    else:
      mockfmodel = fmodel.deep_copy()




    print >> log, "Making new data"
    print >> log, "==============="
    print >> log
    print >> log

    new_data_builder = error_swap( miller_array,
                                   fmodel,
                                   mockfmodel )
    new_data = new_data_builder.new_obs
    # we now have to write the data actually

    print >> log, "Writing new data set"
    print >> log, "===================="

    mtz_dataset = new_data.as_mtz_dataset(
      column_root_label="FOBS")
    mtz_dataset.mtz_object().write(
      file_name=params.simul_utils.output.hklout)




if (__name__ == "__main__" ):
  simul_utils(sys.argv[1:])
