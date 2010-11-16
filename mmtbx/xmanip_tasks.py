from cctbx import miller
import cctbx.xray.structure_factors
from libtbx.utils import Sorry
import iotbx.phil
from iotbx.pdb import xray_structure
from mmtbx.scaling import fa_estimation, pair_analyses, relative_scaling
import sys

master_params = iotbx.phil.parse("""
      task = *get_dano get_diso lsq_scale sfcalc custom None
      .type=choice
      .help="Possible tasks"
      output_label_root=None
      .type=str
      .help="Output label root"
      get_dano
      .help="Get ||F+| - |F-|| from input data."
      {
        input_data = None
        .type=str
      }

      get_diso
      .help="Get |Fder|-|Fnat|"
      {
        native = None
        .type=str
        .help="Name of native data"
        derivative = None
        .type=str
        .help="Name of derivative data"
        use_intensities=True
        .type=bool
        .help="Scale on intensities"
        use_weights=True
        .type=bool
        .help="Use experimental sigmas as weights in scaling"
        scale_weight=True
        .type=bool
        .help="Whether or not to scale the sigmas during scaling"
      }

      lsq_scale{
        input_data_1 = None
        .type=str
        .help="Reference data"
        input_data_2 = None
        .type=str
        .help="Data to be scaled"
        use_intensities=True
        .type=bool
        .help="Scale on intensities"
        use_weights=True
        .type=bool
        .help="Use experimental sigmas as weights in scaling"
        scale_weight=True
        .type=bool
        .help="Whether or not to scale the sigmas during scaling"
      }

      sfcalc{
        fobs = None
        .type=str
        .help = "Data name of observed data"
        output = 2mFo-DFc mFo-DFc *complex_fcalc abs_fcalc intensities
        .type=choice
        .help="Output coefficients"
        use_bulk_and_scale = *as_estimated user_upplied
        .type=choice
        .help = "estimate or use parameters given by user"
        bulk_and_scale_parameters
        .help = "Parameters used in the structure factor calculation. Ignored if experimental data is given"
        {
          d_min = 2.0
          .type=float
          .help = "resolution of the data to be calculated."
          overall
          .help = "Bulk solvent and scaling parameters"
          {
            b_cart
            .help = "Anisotropic B values"
            {
              b_11 = 0
              .type=float
              b_22 = 0
              .type=float
              b_33 = 0
              .type=float
              b_12 = 0
              .type=float
              b_13 = 0
              .type=float
              b_23 = 0
              .type=float
            }
            k_overall = 0.1
            .type=float
            .help = "Overall scalar"
          }
          solvent
          .help = "Solvent parameters"
          {
            k_sol = 0.3
            .type=float
            .help="Solvent scale"
            b_sol = 56.0
            .type=float
            .help="Solvent B"
          }
        }
      }

     custom
     .help = "A custom script that uses miller_array data names as variables."
     {
       code = None
       .help = "A piece of python code"
       .type=str
       show_instructions = True
       .help = "Some instructions"
       .type = bool
     }


      """)


def patch_miller_arrays_as_names( names ):
  result = []
  for name, number in zip(names, range(len(names)) ):
    tmp_result = "%s =  miller_arrays[ %i ].deep_copy()"%(name,number)
    result.append( compile( tmp_result, '<string>', 'exec' )  )

  return result



def get_dano(names, miller_arrays, xray_structure, parameters, out ):
  miller_array = None
  if parameters.input_data is None:
    if len(miller_arrays)==1:
      miller_array = miller_arrays[0]
  else:
    if names.has_key( parameters.input_data ):
      miller_array = miller_arrays[ names[ parameters.input_data ] ]
    else:
      raise Sorry("Unknown data name.")

  if miller_array.is_xray_intensity_array():
    miller_array = miller_array.f_sq_as_f()
  assert miller_array.is_xray_amplitude_array()



  pair_generator = fa_estimation.ano_scaling( miller_array )
  plus  = pair_generator.x1p.deep_copy()
  minus = pair_generator.x1n.deep_copy()
  delta_gen = pair_analyses.delta_generator( plus,
                                             minus )
  deltas = delta_gen.abs_delta_f.deep_copy()
  return deltas

def get_diso(names, miller_arrays, xray_structure, parameters, out):
  #first scale please
  if parameters.native is None:
    raise Sorry("Please define native data name")
  if parameters.derivative is None:
    raise Sorry("Please define derivative data name")

  native=None
  derivative=None

  if names.has_key( parameters.native ):
    native = miller_arrays[ names[parameters.native] ].deep_copy()
  else:
    raise Sorry("Unknown data name: >>%s<<"%(parameters.native) )

  if names.has_key( parameters.derivative ):
    derivative = miller_arrays[ names[parameters.derivative] ].deep_copy()
  else:
    raise Sorry("Unknown data name: >>%s<<"%(parameters.derivative) )

  scaler = relative_scaling.ls_rel_scale_driver(
    miller_native     = native,
    miller_derivative = derivative,
    use_intensities   = parameters.use_intensities,
    scale_weight      = parameters.scale_weight,
    use_weights       = parameters.use_weights)
  #
  scaler.show(out=out)

  if native.is_xray_intensity_array():
    native = native.f_sq_as_f()
  if derivative.is_xray_intensity_array():
    derivative = derivative.f_sq_as_f()

  delta_gen = pair_analyses.delta_generator( derivative,
                                             native )
  deltas = delta_gen.delta_f.deep_copy()
  return deltas

def lsq_scale(names, miller_arrays, xray_structure, parameters, out):
  if parameters.input_data_1 is None:
    raise Sorry("Please define input_data_1")
  if parameters.input_data_2 is None:
    raise Sorry("Please define input_data_2")

  input_data_1 = None
  input_data_2 = None

  if names.has_key( parameters.input_data_1 ):
    input_data_1 = miller_arrays[ names[parameters.input_data_1] ].deep_copy()
  else:
    raise Sorry("Unknown data name: >>%s<<"%(parameters.input_data_1) )

  if names.has_key( parameters.input_data_2 ):
    input_data_2 = miller_arrays[ names[parameters.input_data_2] ].deep_copy()
  else:
    raise Sorry("Unknown data name: >>%s<<"%(parameters.input_data_2) )

  scaler = relative_scaling.ls_rel_scale_driver(
    miller_native     = input_data_1,
    miller_derivative = input_data_2,
    use_intensities   = parameters.use_intensities,
    scale_weight      = parameters.scale_weight,
    use_weights       = parameters.use_weights)
  #
  scaler.show(out=out)
  return scaler.scaled_original_derivative.deep_copy()


def sfcalc(names, miller_arrays, xray_structure, parameters, out):
  from mmtbx import f_model
  f_obs = None
  if parameters.fobs is None:
    if parameters.output not in ["complex_fcalc", "abs_fcalc", "intensities" ]:
      raise Sorry("Experimental data is needed for %s coefficients.\n Please supply Fobs")
    else:
      f_obs = abs(xray_structure.structure_factors(
        d_min          = parameters.bulk_and_scale_parameters.d_min,
        anomalous_flag = False).f_calc())
  else:
    f_obs = miller_arrays[ names[parameters.fobs] ].deep_copy()

  if f_obs.is_xray_intensity_array():
    f_obs = f_obs.f_sq_as_f()

  flags = f_obs.generate_r_free_flags(fraction = 0.1,
                                      max_free = 99999999)
  b_cart = [parameters.bulk_and_scale_parameters.overall.b_cart.b_11,
            parameters.bulk_and_scale_parameters.overall.b_cart.b_22,
            parameters.bulk_and_scale_parameters.overall.b_cart.b_33,
            parameters.bulk_and_scale_parameters.overall.b_cart.b_12,
            parameters.bulk_and_scale_parameters.overall.b_cart.b_13,
            parameters.bulk_and_scale_parameters.overall.b_cart.b_23 ]

  fmodel = f_model.manager( xray_structure   = xray_structure,
                            r_free_flags     = flags,
                            target_name      = "ls_wunit_k1",
                            f_obs            = f_obs,
                            b_cart           = b_cart,
                            k_sol            = parameters.bulk_and_scale_parameters.solvent.k_sol,
                            b_sol            = parameters.bulk_and_scale_parameters.solvent.b_sol)

  if parameters.use_bulk_and_scale == "as_estimated":
    if parameters.fobs is not None:
      fmodel.update_solvent_and_scale(out=out)

  result = None
  if parameters.output in  ["complex_fcalc", "abs_fcalc", "intensities" ]:
    result = fmodel.f_model()
    if parameters.output == "complex_fcalc":
      result = result
    if parameters.output == "abs_fcalc":
      result = abs( result )
    if parameters.output == "intensities":
      result = abs(result).f_as_f_sq()
  else:
    if parameters.output == "2mFo-DFc":
      result = fmodel.electron_density_map().map_coefficients(map_type = "2m*Fobs-D*Fmodel")
    if parameters.output == "mFo-DFc":
      # XXX BUG ?
      result = fmodel.electron_density_map().map_coefficients(map_type = "2m*Fobs-D*Fmodel")

  assert result is not None
  return result



def show_restricted_custom_names(restricted_names, out):
  print >> out, "Restricted data set names are:"
  for rn in restricted_names:
    print >> out, "    -   %s"%(rn)

def print_custom_instructions(out):
  print >> out, "The custom function in the manipulate miller task of xmanip allows one to submit a small (or large)"
  print >> out, "snippet of python code, have it executed and have a single miller array returned and written to file."
  print >> out, "If one is familiar with python and the cctbx in general, this function allows one to quickly perform"
  print >> out, "complex tasks relating reflection files without having the overhead of writing a user interface."
  print >> out, "Data set names given to the miller arrays in the main (user specified) input file, are actual variable names"
  print >> out, "and are stored as a cctbx.miller.array object. A pdb file that was read in, is stored in the object named "
  print >> out, "xray_structure. Note that not many safeguards are in place: make sure your code snippet is proper python!"
  print >> out, "Please note that there are some restriction on variable names: the should not contains spaces or have the name"
  print >> out, "of local variables or functions. By default, a variable named 'result' is returned"



def custom(names, miller_arrays, xray_structure, params, out):

  restricted_names = [ "restricted_names", "names", "miller_arrays", "params", "out",
                       "get_dano", "get_diso", "custom", "sfalc", "patch_miller_arrays_as_names",
                       "lsq_scale", "manipulate_miller", "show_restricted_custom_names", "print_custom_instructions" ]

  if params.show_instructions:
    print_custom_instructions(out)
    show_restricted_custom_names(restricted_names, out)


  #check if all variable names are legal
  for name in names:
    if " " in name:
      raise Sorry("Sorry, no spaces allowed in data set name >%s< to avoid compilation problems."%(name) )
    if name in restricted_names:
      show_restricted_custom_names( restricted_names )
      raise Sorry("The data set name >%s< is restricted to avoid compilation problems." %(name) )

  #first make variables from the names please
  tmp_names = patch_miller_arrays_as_names(names)
  for instruction in tmp_names:
    exec instruction
  result = None
  # now we have to evaulate the code
  print >> out, "Trying to evaluate the code as shown below"
  print >> out, "------------------------------------------"
  print >> out, params.code
  print >> out, "------------------------------------------"
  user_code = compile( params.code, '<string>', 'exec' )
  exec user_code

  return result


def manipulate_miller(names, miller_arrays, xray_structure, params, out=None):
  if out is None:
    out = sys.stdout
  #define a number of function pointers
  function_pointer = {
                       "get_dano" : get_dano,
                       "get_diso" : get_diso,
                       "lsq_scale": lsq_scale,
                       "sfcalc"   : sfcalc,
                       "custom"   : custom,
                     }

  #Now pay attention please
  function_arguments = None
  #these two lines allow me quickly lift the appropriate set of
  #parameters from the file scope without a lengthy set of if statements
  patch = compile("function_arguments = params.%s"%(params.task),'<string>','exec' )
  exec patch
  result = function_pointer[ params.task ]( names,
                                            miller_arrays,
                                            xray_structure,
                                            function_arguments,
                                            out)
  return result
