from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.mpi_integrate
from libtbx.phil import parse

default_steps = [
  'input',
  'balance', # balance input load
  'integrate',
]

from xfel.merging.application.input.file_loader import simple_file_loader
#class noprune_file_loader(simple_file_loader):
def prune_reflection_table_keys(self, reflections):
  return reflections # noop
#import xfel.merging.application.input.file_loader
#xfel.merging.application.input.file_loader.simple_file_loader = noprune_file_loader
simple_file_loader.prune_reflection_table_keys = prune_reflection_table_keys

integrate_phil_str = '''
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope

  dispatch {
    coset = False
      .expert_level = 2
      .type = bool
      .help = Within the integrate dispatcher, integrate a sublattice coset intended to represent \
              negative control spots with no Bragg diffraction.
  }

  integration {
    include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope
    coset {
      transformation = 6
        .type = int(value_min=0, value_max=6)
        .multiple = False
        .help = The index number(s) of the modulus=2 sublattice transformation(s) used to produce distince coset results. \
                0=Double a, 1=Double b, 2=Double c, 3=C-face centering, 4=B-face centering, 5=A-face centering, 6=Body centering \
                See Sauter and Zwart, Acta D (2009) 65:553
    }
    debug {
      delete_shoeboxes = True
        .type = bool
        .help = "Delete shoeboxes immediately before saving files. This option"
                "in combination with debug.output=True enables intermediate"
                "processing steps to make use of shoeboxes."
    }
    integration_only_overrides {
      trusted_range = None
        .type = floats(size=2)
        .help = "Override the panel trusted range (underload and saturation) during integration."
        .short_caption = "Panel trusted range"
    }
  }

  output {
    composite_output = None
      .type = bool
      .expert_level = 4
      .help = Unused
    integrated_filename = None
      .type = str
      .expert_level = 4
      .help = Unused
    integrated_experiments_filename = None
      .type = str
      .expert_level = 4
      .help = Unused
  }
  profile {
    gaussian_rs {
      parameters {
        sigma_b_cutoff = 0.1
          .type = float
          .help = Maximum sigma_b before the image is rejected
      }
    }
  }
'''

program_defaults_phil_str = """
input.keep_imagesets=True
input.parallel_file_load.reset_experiment_id_column = True
"""

from dials.command_line.stills_process import program_defaults_phil_str as dials_program_defaults_phil_str
from xfel.merging.application.phil.phil import dispatch_phil, input_phil, output_phil
phil_scope = parse(integrate_phil_str + dispatch_phil + input_phil + output_phil, process_includes=True).fetch(parse(dials_program_defaults_phil_str)).fetch(parse(program_defaults_phil_str))
import xfel.merging.application.phil.phil
xfel.merging.application.phil.phil.phil_scope = phil_scope

if __name__ == '__main__':
  from xfel.merging.command_line.merge import Script
  script = Script()

  result = script.run()
