from __future__ import absolute_import, division, print_function
import mmtbx.scaling.plan_sad_experiment
from mmtbx.scaling.plan_sad_experiment import get_fp_fdp, get_residues_and_ha
import iotbx.phil
from libtbx.utils import Sorry, null_out
from libtbx import runtime_utils
from libtbx import Auto
import os.path
import sys
from six.moves import range

master_params = """
include scope libtbx.phil.interface.tracking_params
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

  ncs_copies = None
    .type = int
    .short_caption = NCS copies
    .help = Optional estimate of NCS copies in your crystals (only used if \
            a data file is supplied).
   .input_size = 64

  solvent_fraction = None
    .type = float
    .short_caption = Solvent fraction
    .help = Optional estimate of solvent fraction in your crystals (0 to 1)
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

  sites_min = None
    .type = int
    .short_caption = Low bound for sites
    .help = If you set sites_min and sites_max and not sites the sites will\
            be varied from sites_min to sites_max
    .input_size = 64

  sites_max = None
    .type = int
    .short_caption = Upper bound for sites
    .help = If you set sites_min and sites_max and not sites the sites will\
            be varied from sites_min to sites_max
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

   i_over_sigma_range_low = None
     .type = float
     .short_caption = Lower range for i_over_sigma
     .help = If you set i_over_sigma_range_low and i_over_sigma_range_high \
         then the value of i_over_sigma will be varied between these limits

   i_over_sigma_range_high = None
     .type = float
     .short_caption = Upper range for i_over_sigma
     .help = If you set i_over_sigma_range_low and i_over_sigma_range_high \
         then the value of i_over_sigma will be varied between these limits

   steps = 20
     .type = int
     .short_caption = Steps
     .help = Number of steps for sampling ranges (i.e., i_over_sigma_low to \
        i_over_sigma_high)

   min_in_bin = 50
     .type = int
     .short_caption = Minimum point per bin
     .help = Minimum data points per bin in Bayesian estimation. Higher values \
             smooth the predictor.

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

   bayesian_updates = False
     .type = bool
     .short_caption = Bayesian updates
     .help = Use Bayesian updates of half-dataset CC and signal. First \
             predict these values using standard approach, then use empirical \
             half-dataset CC and signal for a training set of datasets to \
             re-estimate these values.  This helps correct for typical errors \
             in measurement and typical resolution-dependent effects.\
             Note that if you use bayesian_updates=True then the predictions \
             may not vary smoothly with resolution or changes in parameters.

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
"""
master_phil = iotbx.phil.parse(master_params, process_includes=True)

def get_params(args,out=sys.stdout):
  command_line = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    reflection_file_def="input_files.data",
    seq_file_def="crystal_info.seq_file")
  params = command_line.work.extract()
  print("\nPlan a SAD experiment\n", file=out)
  master_phil.format(python_object=params).show(out=out)
  return params

def setup_params(params, out):
  if not params.crystal_info.wavelength:
    raise Sorry("Please supply a wavelength for data collection.")
  if not params.crystal_info.atom_type:
    raise Sorry(
      "Please supply an atom_type for the anomalously-scattering atom.")
  if not params.crystal_info.f_double_prime:
     fp_fdp=get_fp_fdp(
      atom_type=params.crystal_info.atom_type,
      wavelength=params.crystal_info.wavelength,out=out)
     if fp_fdp is not None:
       params.crystal_info.f_double_prime=fp_fdp.fdp()
     else:
       raise Sorry(
        "Please specify f_double_prime as the wavelength %7.2f A is" %(
         params.crystal_info.wavelength)+
         "\nout of range of the Sasaki tables used here.")

  if params.crystal_info.solvent_fraction and \
     params.crystal_info.solvent_fraction > 1.01:
    raise Sorry("Solvent fraction should be from 0 to 1")

  if params.crystal_info.seq_file and \
       os.path.isfile(params.crystal_info.seq_file):
    residues,sites,number_of_s,solvent_fraction,ncs_copies=get_residues_and_ha(
      seq_file=params.crystal_info.seq_file,
      atom_type=params.crystal_info.atom_type,
      chain_type=params.crystal_info.chain_type,
      solvent_fraction=params.crystal_info.solvent_fraction,
      data=params.input_files.data,
      ncs_copies=params.crystal_info.ncs_copies,
      out=out)
    if not params.crystal_info.residues:
      print("Number of residues based on sequence file: %d" %(
        residues), file=out)
      params.crystal_info.residues=residues
    if not params.crystal_info.number_of_s:
      print("Number of S atoms based on sequence file: %d" %(
        number_of_s), file=out)
      params.crystal_info.number_of_s=number_of_s

    if not params.crystal_info.sites:
      print("Number of sites for anomalously-scattering atom "+\
        "based on sequence file: %d" %( sites), file=out)
      params.crystal_info.sites=sites

    if ncs_copies and not params.crystal_info.ncs_copies:
      print("NCS copies "+\
        "based on sequence file and data : %d" %( ncs_copies), file=out)
      params.crystal_info.ncs_copies=ncs_copies

    if solvent_fraction and not params.crystal_info.solvent_fraction:
      print("Solvent fraction "+\
        "based on sequence file and data : %5.2f" %( solvent_fraction), file=out)
      params.crystal_info.solvent_fraction=solvent_fraction

  else:
    if params.crystal_info.number_of_s is None and \
         params.include_weak_anomalous_scattering is True:
      raise Sorry("Sorry need a sequence file or number_of_s if "+
        "\ninclude_weak_anomalous_scattering=True")
    elif params.crystal_info.number_of_s is None and \
         params.include_weak_anomalous_scattering is Auto:
      print("Note: not applying include_weak_anomalous_scattering as"+\
        " no sequence \nfile or number_of_s are supplied", file=out)
      params.include_weak_anomalous_scattering=False

  if params.crystal_info.solvent_fraction is None:
    params.crystal_info.solvent_fraction=0.50  # just guess

  if not params.crystal_info.residues:
    raise Sorry("Please specify number of residues (residues=500) or a sequence file")
  if params.crystal_info.sites:
    pass # OK
  elif (params.crystal_info.sites_min and params.crystal_info.sites_max):
    pass # OK
  else:
    raise Sorry(
      "Please specify number of sites or a sequence file and atom_type")


class result_table:
  def __init__(self):
    self.table_rows=[]
    self.table_header=[]
    self.number_of_columns=0

  def add_table_header(self,header): # must be first
    self.table_header=header
    self.number_of_columns=len(header)

  def add_table_row(self,row):
    assert self.table_header and len(row)==len(self.table_header)
    self.table_rows.append(row)

  def get_formats(self,buffer=0):
    self.widths=[]
    self.formats=[]
    for i in range(self.number_of_columns):
      w=len(self.table_header[i])
      for tr in self.table_rows:
        w=max(w,len(tr[i]))
      self.widths.append(w)
      self.formats.append("%s%ss" %("%",w+buffer))

  def show_summary(self,buffer=3,gui_output=False,out=sys.stdout):
    assert not gui_output # not implemented yet
    self.get_formats(buffer=buffer)
    for i in range(self.number_of_columns):
      print(self.formats[i] %(self.table_header[i]), end=' ', file=out)
    print(file=out)

    for tr in self.table_rows:
      for i in range(self.number_of_columns):
        print(self.formats[i] %(tr[i]), end=' ', file=out)
      print(file=out)

def run(args,params=None,return_plan=False,out=sys.stdout):
  # NOTE: can call with params and skip reading any files.
  if not params:
    params=get_params(args,out=out)

  setup_params(params, out=out)

  if params.crystal_info.sites_min and params.crystal_info.sites_max:
    return run_varying_sites(params,out=out)
  elif params.i_over_sigma_range_low and params.i_over_sigma_range_high:
    return run_varying_i_over_sigma(params,out=out)
  else:
    return run_with_params(params,out=out)

def run_varying_i_over_sigma(params,out=sys.stdout):
  from copy import deepcopy
  local_params=deepcopy(params)
  local_params.i_over_sigma_range_low=None
  local_params.i_over_sigma_range_high=None
  local_params.control.fixed_resolution=True

  local_params.control.show_summary=True
  t=result_table()
  t.add_table_header([
     "I/sigI",
     "cc_half",
     "cc*_anom",
     "Signal",
     "p(Substr)",
     "FOM",
   ])

  delta=max(0.001,
    (params.i_over_sigma_range_high-params.i_over_sigma_range_low)/max(1,
     params.steps))
  i_over_sigma=params.i_over_sigma_range_low
  while i_over_sigma <= params.i_over_sigma_range_high+0.01:
    local_params.i_over_sigma=i_over_sigma
    plan=run_with_params(local_params,quiet=True,out=out)
    [dmin,nsites,nrefl,fpp,local_i_over_sigma,
        sigf,cc_half_weak,cc_half,cc_ano_weak,cc_ano,s_ano,solved,fom]=\
      plan.representative_values
    t.add_table_row([
        "%6.2f" % (i_over_sigma),
        "%6.3f" % (cc_half),
        "%6.3f" % ( cc_ano),
        "%5.2f" % ( s_ano),
        "%5.2f" % (solved),
        "%4.3f" % (fom),
      ])
    i_over_sigma+=delta

  plan.show_characteristics(out=out)
  print("\nExpected data utility varying the value of overall I/sigI", file=out)
  t.show_summary(out=out)

def run_varying_sites(params,out=sys.stdout):
  if not params.i_over_sigma:
    raise Sorry("For varying sites you need to set i_over_sigma")
  from copy import deepcopy
  local_params=deepcopy(params)
  local_params.crystal_info.sites_min=None
  local_params.crystal_info.sites_max=None
  local_params.control.show_summary=True
  local_params.control.fixed_resolution=True
  t=result_table()
  t.add_table_header([
     "Sites",
     "cc_half",
     "cc*_anom",
     "Signal",
     "p(Substr)",
     "FOM",
   ])

  if params.crystal_info.sites_min>params.crystal_info.sites_max:
    raise Sorry("Please set sites_min < sites_max")
  for sites in range(params.crystal_info.sites_min,
     params.crystal_info.sites_max+1):
    local_params.crystal_info.sites=sites
    plan=run_with_params(local_params,quiet=True,out=out)
    [dmin,nsites,nrefl,fpp,local_i_over_sigma,
        sigf,cc_half_weak,cc_half,cc_ano_weak,cc_ano,s_ano,solved,fom]=\
      plan.representative_values
    t.add_table_row([
        "%3d" % (nsites),
        "%6.3f" % (cc_half),
        "%6.3f" % ( cc_ano),
        "%5.2f" % ( s_ano),
        "%5.1f" % (solved),
        "%4.3f" % (fom),
      ])
  plan.show_characteristics(out=out)

  print("\nExpected data utility varying the number of sites", file=out)
  t.show_summary(out=out)


def run_with_params(params,quiet=False,out=sys.stdout):
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
    min_in_bin=params.min_in_bin,
    data=params.input_files.data,
    data_labels=params.input_files.data_labels,
    resolution=params.crystal_info.resolution,
    b_value=params.crystal_info.b_value,
    b_value_anomalous=params.crystal_info.b_value_anomalous,
    fixed_resolution=params.control.fixed_resolution,
    occupancy=params.crystal_info.occupancy,
    ideal_cc_anom=params.ideal_cc_anom,
    bayesian_updates=params.bayesian_updates,
    include_weak_anomalous_scattering=params.include_weak_anomalous_scattering,
    intrinsic_scatterers_as_noise=params.intrinsic_scatterers_as_noise,)
  if quiet:
    return plan
  elif params.control.show_summary:
    plan.show_summary()
  else:
    return plan.show(out=out)

def validate_params(params):
  return setup_params(params, out=null_out())

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    return run(self.args, out=sys.stdout)

def finish_job(result):
  return ([], [])

if __name__=="__main__":
  run(sys.argv[1:])

