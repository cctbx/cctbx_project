
from __future__ import division
import mmtbx.scaling.plan_sad_experiment
from mmtbx.scaling.plan_sad_experiment import get_fpp, get_residues_and_ha
import iotbx.phil
from libtbx.utils import Sorry, null_out
from libtbx import runtime_utils
import os.path
import sys

master_phil = iotbx.phil.parse("""
crystal_info {
  resolution = None
    .type = float
    .help = Optional high-resolution limit. If specified, the calculation \
            is only carried out at this resolution
    .short_caption = High-resolution limit
    .style = resolution
   .input_size = 64

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
}
target_signal = 30.
  .type = float
  .short_caption = Target anomalous signal
  .help = The anomalous signal that you would like to obtain. \
          The value of I/sigma will be adjusted to obtain this signal.\
          Typically you will need a signal of 15-30 so solve the \
          substructure.
  .input_size = 64
min_cc_ano = 0.15
  .type = float
  .short_caption = Target minimum anomalous correlation
  .help = You can set the target minimum (true) anomalous correlation \
          (CC*_ano). This value affects the phasing accuracy after \
          the substructure is determined.
  .input_size = 64

control {
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
     params.crystal_info.f_double_prime=get_fpp(
      atom_type=params.crystal_info.atom_type,
      wavelength=params.crystal_info.wavelength,out=out)

  if params.crystal_info.seq_file and \
       os.path.isfile(params.crystal_info.seq_file):
    residues,sites=get_residues_and_ha(params.crystal_info.seq_file,
      params.crystal_info.atom_type,
      params.crystal_info.chain_type,out=out)
    if not params.crystal_info.residues:
      print >>out,"Number of residues based on sequence file: %d" %(
        residues)
      params.crystal_info.residues=residues

    if not params.crystal_info.sites:
      print >>out,"Number of sites for anomalously-scattering atom "+\
        "based on sequence file: %d" %( sites)
      params.crystal_info.sites=sites

  if not params.crystal_info.residues:
    raise Sorry("Please specify number of residues or a sequence file")
  if not params.crystal_info.sites:
    raise Sorry(
      "Please specify number of sites or a sequence file and atom_type")

def run(args,params=None,out=sys.stdout):
  # NOTE: can call with params and skip reading any files.
  if not params:
    params=get_params(args,out=out)

  setup_params(params, out=out)

  return mmtbx.scaling.plan_sad_experiment.estimate_necessary_i_sigi(
    chain_type=params.crystal_info.chain_type,
    residues=params.crystal_info.residues,
    solvent_fraction=params.crystal_info.solvent_fraction,
    nsites=params.crystal_info.sites,
    fpp=params.crystal_info.f_double_prime,
    target_s_ano=params.target_signal,
    min_cc_ano=params.min_cc_ano,
    dmin=params.crystal_info.resolution).show(out=out)

def validate_params (params) :
  return setup_params(params, out=null_out())

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(self.args, out=sys.stdout)

def finish_job (result) :
  return ([], [])

if __name__=="__main__":
  run(sys.argv[1:])
