# LIBTBX_SET_DISPATCHER_NAME phenix.plan_sad_experiment

from __future__ import division
import mmtbx.scaling.plan_sad_experiment
from mmtbx.scaling.plan_sad_experiment import get_fp_fdp, get_residues_and_ha
import iotbx.phil
from libtbx.utils import Sorry, null_out
from libtbx import runtime_utils
from libtbx import Auto
import os.path
import sys

master_phil = iotbx.phil.parse("""
input_files {
    data = None
      .type = path
      .help = Data file (I or I+ and I- or F or F+ and F-).  \
         Any standard format is fine. 
      .short_caption = Data file 
      .style = bold file_type:hkl input_file process_hkl child:fobs:data_labels\
        child:space_group:space_group child:unit_cell:unit_cell anom
    data_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of anomalous data to use.\
        Not necessary if your input file has only one set of anomalous data.
      .short_caption = Data label
      .style = bold renderer:draw_fobs_label_widget
}
crystal_info {
  resolution = None
    .type = float
    .help = High-resolution limit.  \
       Either a high-resolution limit or a \
       (Wilson) b_value or both is required
    .short_caption = High-resolution limit
    .style = resolution
   .input_size = 64

  b_value = None
    .type = float
    .help = Estimated Wilson B-value for the dataset. \
       Either a high-resolution limit or a \
       (Wilson) b_value or both is required

  b_value_anomalous = None
    .type = float
    .help = Estimated Wilson B-value for the anomalously-scattering atoms. \
       Normally leave as None and it will be estimated from b_value.

  seq_file = None
    .type = path
    .help = "Optional sequence file (1-letter code)."
             "Separate chains with a blank line or line starting with &gt;."
    .short_caption = Sequence file
    .style = bold file_type:seq input_file

  chain_type = *PROTEIN RNA DNA
    .type = choice
    .short_caption = Chain type
    .help = Chain type (PROTEIN RNA DNA). This is used to estimate the \
            number of atoms from the number of residues

  solvent_fraction = 0.5
    .type = float
    .short_caption = Solvent fraction
    .help = Optional estimate of solvent fraction in your crystals
   .input_size = 64

  residues = None
    .type = int
    .short_caption = Number of residues
    .help = The number of residues in the molecule or asymmetric unit. \
          Note that it is the ratio of residues to anomalously-scattering \
          atoms that matters.
    .input_size = 64

  atom_type = None
    .type = str
    .short_caption = Anomalously-scattering atom
    .help = Optional name of anomalously-scattering atom.  If supplied, \
            you also need to supply the wavelength for X-ray data \
            collection.  If not supplied, then you need to supply a \
            value for f_double_prime.
    .style = renderer:draw_phaser_scatterer_widget
    .input_size = 64

  number_of_s = None
    .type = int
    .short_caption = Number of S atoms
    .help = You can specify the number of S atoms in the asymmetric unit. \
            Only used if include_weak_anomalous_scattering=True.  If not \
            set, the number is guessed from the sequence file if present.
  f_double_prime = None
    .type = float
    .short_caption = F-double-prime
    .help = F-double-prime value for the anomalously-scattering atom. \
            Alternatively you can specify the atom type and wavelength.
    .input_size = 64

  wavelength = None
    .type = float
    .short_caption = Wavelength
    .help = Wavelength for X-ray data collection.  If supplied, also \
            specify the atom_type. \
            Alternatively you can specify the value of f_double_prime
    .input_size = 64

  sites = None
    .type = int
    .short_caption = Number of anomalously-scattering atoms
    .help = The number of anomalously-scattering atoms in the molecule \
          or asymmetric unit. \
          Note that it is the ratio of residues to anomalously-scattering \
          atoms that matters.
    .input_size = 64

  occupancy = 1
    .type = float
    .short_caption = Occupancy of anomalously-scattering atoms
    .help = Estimate of occupancy of anomalously-scattering atoms


   }
   i_over_sigma = None
     .type = float
     .short_caption = I/sigI
     .help = Optional I/sigI.  If supplied, \
        the expected values of half-dataset correlation and cc*_ano based \
        on this I/sigI be calculated.

   max_i_over_sigma = 100
     .type = float
     .short_caption = Maximum I/sigI
     .help = Limit search of necessary I/sigI to less than this value.  \
             You might increase this if you plan to do a very careful or very \
             high-multiplicity experiment.

   target_signal = 30.
       .type = float
       .short_caption = Target anomalous signal
       .help = The anomalous signal that you would like to obtain. \
               The value of I/sigma will be adjusted to obtain this signal.\
               Typically you will need a signal of 15-30 so solve the \
               substructure.
   min_cc_ano = 0.15
     .type = float
     .short_caption = Target minimum anomalous correlation
     .help = You can set the target minimum (true) anomalous correlation \
             (CC*_ano). This value affects the phasing accuracy after \
             the substructure is determined.

   ideal_cc_anom = 0.75
     .type = float
     .short_caption = Anomalous correlation with perfect data
     .help = The ideal_cc_anom is the expected anomalous \
             correlation between an accurate model with isotropic anomalous \
             scatterers and perfectly-measured data. The ideal_cc_anom is  \
             determined empirically.  It is typically not unity because \
             anomalous scatterers may have multiple locations with low \
             occupancy or may be non-isotropic. A value of about 0.75 \
             is a reasonable guess.

   include_weak_anomalous_scattering = Auto
     .type = bool
     .short_caption = Include weak anomalous scattering
     .help = At longer wavelengths the scattering of C, N, and O become \
             significant relative to S. Default is to consider the \
             scattering from C, N, O as noise.  Additionally, \
             (see intrinsic_scatterers_as_noise) if intrinsic \
             anomalous scatterers (P and S) are weak, they will be counted \
             as noise.  This weak anomalous scattering is effectively noise \
             and has the same effect as the ideal_cc_anom but it can be \
             calculated from the composition. Its effects are added to those \
             modeled by the ideal_cc_anom parameter. Default is to include \
             weak anomalous scattering if a sequence file or the number of \
             sulfurs is provided

   intrinsic_scatterers_as_noise = None
     .type = bool
     .short_caption = Intrinsic scatterers as noise
     .help = Applies if include_weak_anomalous_scattering=True.\
             You can choose to treat any intrinsic scatterers (S for \
             protein, P for nucleic acid) as noise, just like any \
             contributions from C, N, or O atoms. This is default if \
             anomalous scattering (f-double-prime) from these atoms is \
             less than half that of your specified anomalous scatterer. \
             Alternatively these atoms are excluded from the noise \
             calculation and are assumed to be included in the \
             number of sites you specify.

   bayesian_estimates = True
     .type = bool
     .short_caption = Bayesian estimates
     .help = Use Bayesian estimates of half-dataset CC and signal. First \
             predict these values using standard approach, then use empirical \
             half-dataset CC and signal for a training set of datasets to \
             re-estimate these values.  This helps correct for typical errors \
             in measurement and typical resolution resolution-dependent effects.

   control {
      fixed_resolution = False
        .type = bool
        .help = Only run calculation at high_resolution limit
        .short_caption = Run at high_resolution only

      show_summary = False
        .type = bool
        .help = Show summary only
        .short_caption = Show summary only

      verbose = False
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output
   }
""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):
  command_line = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    seq_file_def="crystal_info.seq_file")
  params = command_line.work.extract()
  print >>out,"\nPlan a SAD experiment\n"
  master_phil.format(python_object=params).show(out=out)
  return params

def setup_params (params, out) :
  if not params.crystal_info.f_double_prime:
     params.crystal_info.f_double_prime=get_fp_fdp(
      atom_type=params.crystal_info.atom_type,
      wavelength=params.crystal_info.wavelength,out=out).fdp()

  if params.crystal_info.seq_file and \
       os.path.isfile(params.crystal_info.seq_file):
    residues,sites,number_of_s=get_residues_and_ha(
      params.crystal_info.seq_file,
      params.crystal_info.atom_type,
      params.crystal_info.chain_type,out=out)
    if not params.crystal_info.residues:
      print >>out,"Number of residues based on sequence file: %d" %(
        residues)
      params.crystal_info.residues=residues
    if not params.crystal_info.number_of_s:
      print >>out,"Number of S atoms based on sequence file: %d" %(
        number_of_s)
      params.crystal_info.number_of_s=number_of_s

    if not params.crystal_info.sites:
      print >>out,"Number of sites for anomalously-scattering atom "+\
        "based on sequence file: %d" %( sites)
      params.crystal_info.sites=sites
  else:
    if params.crystal_info.number_of_s is None and \
         params.include_weak_anomalous_scattering is True:
      raise Sorry("Sorry need a sequence file or number_of_s if "+
        "\ninclude_weak_anomalous_scattering=True")
    elif params.crystal_info.number_of_s is None and \
         params.include_weak_anomalous_scattering is Auto:
      print >>out,"Note: not applying include_weak_anomalous_scattering as"+\
        " no sequence \nfile or number_of_s are supplied"
      params.include_weak_anomalous_scattering=False

  if not params.crystal_info.residues:
    raise Sorry("Please specify number of residues or a sequence file")
  if not params.crystal_info.sites:
    raise Sorry(
      "Please specify number of sites or a sequence file and atom_type")

def run(args,params=None,return_plan=False,out=sys.stdout):
  # NOTE: can call with params and skip reading any files.
  if not params:
    params=get_params(args,out=out)

  setup_params(params, out=out)

  plan=mmtbx.scaling.plan_sad_experiment.estimate_necessary_i_sigi(
    chain_type=params.crystal_info.chain_type,
    residues=params.crystal_info.residues,
    number_of_s=params.crystal_info.number_of_s,
    solvent_fraction=params.crystal_info.solvent_fraction,
    nsites=params.crystal_info.sites,
    wavelength=params.crystal_info.wavelength,
    atom_type=params.crystal_info.atom_type,
    fpp=params.crystal_info.f_double_prime,
    target_s_ano=params.target_signal,
    i_over_sigma=params.i_over_sigma,
    max_i_over_sigma=params.max_i_over_sigma,
    min_cc_ano=params.min_cc_ano,
    data=params.input_files.data,
    data_labels=params.input_files.data_labels,
    resolution=params.crystal_info.resolution,
    b_value=params.crystal_info.b_value,
    b_value_anomalous=params.crystal_info.b_value_anomalous,
    fixed_resolution=params.control.fixed_resolution,
    occupancy=params.crystal_info.occupancy,
    ideal_cc_anom=params.ideal_cc_anom,
    bayesian_estimates=params.bayesian_estimates,
    include_weak_anomalous_scattering=params.include_weak_anomalous_scattering,
    intrinsic_scatterers_as_noise=params.intrinsic_scatterers_as_noise,)
  if params.control.show_summary:
    plan.show_summary()
  else:
    return plan.show(out=out)

def validate_params (params) :
  return setup_params(params, out=null_out())

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(self.args, out=sys.stdout)

def finish_job (result) :
  return ([], [])

if __name__=="__main__":
  run(sys.argv[1:])
