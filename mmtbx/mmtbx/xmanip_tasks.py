from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
import cctbx.xray.structure_factors
from cctbx.eltbx.xray_scattering import wk1995
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.math import matrix
from cctbx import adptbx
from mmtbx.scaling import fa_estimation, pair_analyses, relative_scaling
import sys, os

master_params = iotbx.phil.parse("""
      task = *get_dano get_diso lsq_scale None
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

      """)

def get_dano(names, miller_arrays, parameters, out ):
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

def get_diso(names, miller_arrays, parameters, out):
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

def lsq_scale(names, miller_arrays, parameters, out):
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



def manipulate_miller(names, miller_arrays, params, out=None):
  if out is None:
    out = sys.stdout
  #define a number of function pointers
  function_pointer = { "get_dano": get_dano,
                       "get_diso": get_diso,
                       "lsq_scale": lsq_scale }
  #Now pay attention please
  function_arguments = None
  #these two lines allow me quickly lift the appropriate set of
  # parameters from the file scope without a length y set o f if statements
  patch = compile("function_arguments = params.%s"%(params.task),'<string>','exec' )
  exec patch

  result = function_pointer[ params.task ]( names,
                                            miller_arrays,
                                            function_arguments,
                                            out)
  return result
