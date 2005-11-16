import os,sys, string, iotbx.phil

class phil_lego(object):
  """
This class facilitates the construction of phil parameter files
for the FA estimation program FATSO.
"""
  def __init__(self):
    self.scaling_input = """ scaling.input{
__REPLACE__

expert_level=0
.type=int
.expert_level=1
}
"""
    self.basic_info = """basic{
  n_residues=None
  .type=float
  n_bases=None
  .type=float
  n_copies_per_asu=None
  .type=float
}
"""
    self.xray_data_basic="""xray_data{
  unit_cell=None
  .type=unit_cell

  space_group=None
  .type=space_group

  __REPLACE__
}
"""

    self.data_type="""__REPLACE__{
  file_name=None
  .type=path
  labels=None
  .type=strings
}
"""
    self.scaling_strategy="""scaling_strategy
.expert_level=__EXPERT_LEVEL__
{
  __REPLACE__
}
"""

    self.pre_scaler_protocol="""pre_scaler_protocol
.expert_level=__EXPERT_LEVEL__
{
high_resolution=None
.type=float
low_resolution=None
.type=float
aniso_correction=True
.type=bool
outlier_level_wilson=1e-6
.type=float
 outlier_level_extreme=1e-2
.type=float
}"""

    self.scale_protocol="""__REPLACE__
.expert_level=__EXPERT_LEVEL__
{
           target = ls loc *ls_and_loc None
         .type=choice
         iterations = *auto specified_by_max_iterations
         .type=choice
         max_iterations = 2
         .type=int

         least_squares_options{
           use_experimental_sigmas=True
           .type=bool
           scale_data=*intensities amplitudes
           .type=choice
           scale_target=basic *fancy
           .type=choice
         }

         local_scaling_options{
           use_experimental_sigmas=True
           .type=bool
           scale_data=*intensities amplitudes
           .type=choice
           scale_target=*local_moment local_lsq
           .type=choice
           max_depth=10
           .type=int
           target_neighbours=100
           .type=int
           neighbourhood_sphere=1
           .type=int
         }

         outlier_rejection_options{
           cut_level_sigma=3
           .type=float
           cut_level_rms_primary=4
           .type=float
           cut_level_rms_secondary=4
           .type=float
           protocol=solve rms *rms_and_sigma
           .type=choice
         }


}"""

    self.output="""output
{
     log = 'fatso.log'
     .type = path
     hklout = 'fatso.mtz'
     .type = path
     outlabel = '_ATSO'
     .type = str

}
"""


  def default_sad(self):
    outer_level = self.scaling_input
    basic = self.basic_info
    data = self.data_type.replace( '__REPLACE__',
                                      'reference' )
    data = self.xray_data_basic.replace('__REPLACE__',
                                           data )
    scaler = self.pre_scaler_protocol + \
             self.scale_protocol.replace('__REPLACE__',
                                         'ano_protocol' )
    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )
    scaler = self.scaling_strategy.replace('__REPLACE__',
                                           scaler )
    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )
    scaler = scaler.replace( 'ls loc *ls_and_loc None',
                              '*loc None' )
    output = self.output

    result = outer_level.replace('__REPLACE__',
                                 basic+data+scaler+output)
    return result

  def default_sir(self):
    outer_level = self.scaling_input
    basic = self.basic_info
    data = self.data_type.replace( '__REPLACE__',
                                      'native' ) \
                                      + \
            self.data_type.replace( '__REPLACE__',
                                      'derivative' )


    data = self.xray_data_basic.replace('__REPLACE__',
                                           data )


    scaler = self.pre_scaler_protocol + \
             self.scale_protocol.replace('__REPLACE__','iso_protocol' )

    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )
    scaler = self.scaling_strategy.replace('__REPLACE__',
                                           scaler )
    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )

    output = self.output

    result = outer_level.replace('__REPLACE__',
                                 basic+data+scaler+output)
    return result



  def default_siras(self):
    outer_level = self.scaling_input
    basic = self.basic_info
    data = self.data_type.replace( '__REPLACE__',
                                      'native' ) \
                                      + \
            self.data_type.replace( '__REPLACE__',
                                      'derivative' )

    data = self.xray_data_basic.replace('__REPLACE__',
                                           data )

    scaler = self.scale_protocol.replace('__REPLACE__',
                                         'ano_protocol' )
    scaler = scaler.replace('ls loc *ls_and_loc None',
                            '*loc None' )

    scaler = self.pre_scaler_protocol + scaler + \
             self.scale_protocol.replace('__REPLACE__','iso_protocol' )

    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )
    scaler = self.scaling_strategy.replace('__REPLACE__',
                                           scaler )
    scaler = scaler.replace('__EXPERT_LEVEL__',
                            '1' )
    output = self.output

    result = outer_level.replace('__REPLACE__',
                                 basic+data+scaler+output)
    return result




def run(args):
  tester = phil_lego()

  print " ---------- SAD ----------"
  master_params = iotbx.phil.parse( tester.default_sad() )
  master_params.show(expert_level = int(args[0]) )

  print " ---------- SIR ----------"
  master_params = iotbx.phil.parse( tester.default_sir() )
  master_params.show(expert_level=int(args[0]))
  print " ---------- SIRAS ----------"
  master_params = iotbx.phil.parse( tester.default_siras() )
  master_params.show(expert_level=int(args[0]))



if (__name__ == "__main__"):
  run(sys.argv[1:])
