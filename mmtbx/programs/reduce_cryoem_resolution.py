from __future__ import absolute_import, division, print_function
import os
import math

# Import necessary components for the Program Template and data handling
try:
    from phenix.program_template import ProgramTemplate
    inside_phenix = True
except ImportError:
    # Fallback if running outside a full Phenix environment, assumes libtbx is available
    inside_phenix = False
    from libtbx.program_template import ProgramTemplate

# import cctbx components
from libtbx.utils import Sorry
from libtbx import group_args
from iotbx.cli_parser import run_program
from scitbx.array_family import flex

# import numpy and minimize from scipy
import numpy as np
from scipy.optimize import minimize

# -----------------------------------------------------------------------------
# 1. Define PHIL Parameters (master_phil_str)
# -----------------------------------------------------------------------------

master_phil_str = """
map_model {
  include scope iotbx.map_model_manager.map_phil_str
}

d_min = None
  .type = float
  .help = Nominal resolution of input map(s) (Angstroms). Optional to see initial FSC
  .short_caption = Nominal resolution of input map(s)
  .optional = True
target_d_min = None
  .type = float
  .help = Desired high-resolution limit (Angstroms)
  .optional = False
mute = False
  .type = bool
  .help = Flag to mute standard output
  .optional = True
output_files {
  output_dir = "."
    .type = path
    .help = Directory where output files will be written
    .optional = True
  file_root = "reduced_resolution"
    .type = path
    .help = Filename root for resulting output map file(s)
    .optional = True
  plot = False
    .type = bool
    .help = Flag to make and save FSC plot
    .optional = True
}
"""

# -----------------------------------------------------------------------------
# 2. Define required functions
# -----------------------------------------------------------------------------

def apply_spherical_mask(mmm, d_min):
  unit_cell = mmm.map_manager().unit_cell()
  ucpars = unit_cell.parameters()
  mask_radius = min(ucpars[0], ucpars[1], ucpars[2])/4.
  working_mmm = mmm.deep_copy()
  working_mmm.create_spherical_mask(soft_mask_radius=d_min,
                            mask_radius = mask_radius+d_min)
  working_mmm.apply_mask_to_maps()
  return working_mmm

def FSCfunction(signalratio, deltaB, d_inv):
  sbterm = signalratio * math.exp(-deltaB*d_inv**2/4)
  FSCcalc = sbterm / (1. + sbterm)
  return FSCcalc

def dFSCsqr(x, fsc, d_inv):
  if len(x) != 2:
    raise Sorry("Must have two parameters in dFSCsqr")
  (signalratio, deltaB) = x
  ndat = len(fsc)
  if len(d_inv) != ndat:
    raise Sorry("Mismatch in fsc and d_inv data")
  sumsqr = 0.
  for i in range(ndat):
    FSCcalc = FSCfunction(signalratio, deltaB, d_inv[i])
    sumsqr += (fsc[i] - FSCcalc)**2
  return sumsqr

def fsc_params_from_d_min(d_min):
  signalratio = 34.4680 * math.exp(3.00533/d_min + 4.27895/d_min**2)
  deltaB = 17.1158 + 12.0213*d_min + 21.3225*d_min**2
  return group_args(signalratio=signalratio,
                    deltaB=deltaB)

def d_min_from_fsc_analytical(mc1, mc2, guess_d_min = 3.0):
  fsc_results = mc1.fsc(other=mc2)
  fsc_bins = fsc_results.fsc
  d_inv = fsc_results.d_inv
  start_params = fsc_params_from_d_min(guess_d_min)
  start_sr = start_params.signalratio
  start_dB = start_params.deltaB
  x0 = np.array((start_sr, start_dB))
  initial_simplex = np.array((x0, x0+np.array((start_sr/10.,0.)), x0+np.array((0.,start_dB/10.))))
  results_min = minimize(dFSCsqr, x0, method='Nelder-Mead',
                        args = (fsc_bins,d_inv),
                        bounds = ((1,1000),
                                  (1,10000)),
                        options={'initial_simplex': initial_simplex})
  (signalratio, deltaB) = results_min.x
  if deltaB > 0.:
    d_min = math.sqrt(deltaB/(1.790593 + math.log(signalratio))) / 2.
  else:
    d_min = guess_d_min
  return group_args(d_min=d_min,
                    d_inv=d_inv,
                    signalratio=signalratio,
                    deltaB=deltaB,
                    fsc_bins=fsc_bins)

# -----------------------------------------------------------------------------
# 3. Program Definition
# -----------------------------------------------------------------------------

class Program(ProgramTemplate):
    description = """
    Takes two cryo-EM half-maps directly OR copied from a single full map,
    then adds appropriate noise to half-maps to limit the resulting maps
    to a specified target resolution. Filters output maps to that resolution.
    Initial estimate of d_min is optional but needed to deduce resolution from starting FSC.
    If the plot keyword is set to True, a plot of FSC curves will be generated and saved.

    Usage Example (Half-maps):
      python reduce_cryoEM_resolution.py hm1.mrc hm2.mrc d_min=3.0 target_d_min=4.0 plot=True
    Usage Example (Full map): python reduce_cryoEM_resolution.py full.mrc target_d_min=4.0
    """

    # Required attributes for Program class
    master_phil_str = master_phil_str
    datatypes = ['real_map', 'phil']

    def validate(self):
        """
        Check for required inputs: valid map configuration (2 half-maps OR 1 full map)
        and both resolution parameters.
        """
        # 1. Check map input exclusivity and validity
        nmaps = len(self.data_manager.get_real_map_names())
        if nmaps not in [1,2]:
          raise Sorry("Either one full map or two half maps should be provided")
        if nmaps == 1:
          self.halfmaps = False
        else:
          self.halfmaps = True

        # Make map_model_manager from input map(s)
        self.mmm = self.data_manager.get_map_model_manager(from_phil=True)

        p = self.params
        # 2. Check for required target resolution parameter
        if not p.target_d_min:
            # The target resolution acts as the high resolution limit (d_min) for the filter.
            raise Sorry("Must supply the target output resolution limit (target_d_min).")
        if p.d_min is not None:
          if p.target_d_min < p.d_min:
            raise Sorry("target_d_min should be lower than d_min")

        return True

    def print_version_info(self):
      # Print version info
      import time
      if inside_phenix:
        print ("\n"+60*"*"+"\n"+"  PHENIX reduce_cryoEM_resolution"+
          "  "+str(time.asctime())+"\n"+60*"*"+"\n", file = self.logger)
      else:
        print ("\n"+60*"*"+"\n"+"    reduce_cryoEM_resolution"+
          "  "+str(time.asctime())+"\n"+60*"*"+"\n", file = self.logger)
      print ("Working directory: ", os.getcwd(), "\n", file = self.logger)
      if inside_phenix:
        print ("PHENIX VERSION: ", os.environ.get('PHENIX_VERSION', 'svn'), "\n",
        file = self.logger)

    def run(self):
        """
        Main functionality: Load 1 or 2 maps, set nominal resolution,
        add half-map noise and filter to target resolution, write the output maps.
        """

        params = self.params
        mute = params.mute
        if not mute:
          # print version and date
          self.print_version_info()
        plot = params.output_files.plot
        file_root = params.output_files.file_root
        target_d_min = params.target_d_min
        d_min = params.d_min
        if d_min is None:
          d_min = target_d_min * 0.8
          d_min_given = False
        else:
          d_min_given = True
        log = self.logger
        dm = self.data_manager
        halfmaps = self.halfmaps
        mmm = self.mmm
        if not halfmaps:
          mmm.add_map_manager_by_id(mmm.map_manager(),'map_manager_1')
          mmm.add_map_manager_by_id(mmm.map_manager(),'map_manager_2')

        mc1 = mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_1')
        working_mmm = apply_spherical_mask(mmm, target_d_min)
        mask_fraction = working_mmm.mask_info().mean # Use to correct noise for full box
        wmc1 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_1')
        zeros_array = flex.complex_double(mc1.size(), 0.+0.j)
        mc1_new = mc1.customized_copy(data=zeros_array)
        mc2_new = mc1_new.deep_copy()
        if halfmaps:
          mc2 = mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_2')
          wmc2 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_2')
          if d_min_given:
            d_min_results = d_min_from_fsc_analytical(wmc1,wmc2,guess_d_min=d_min)
            d_min_from_fsc = d_min_results.d_min
            start_d_inv = d_min_results.d_inv
            start_FSC = d_min_results.fsc_bins
            if not mute:
              print("Initial estimate of d_min: ", d_min_from_fsc, file=log)
            if d_min_from_fsc >= target_d_min:
              raise Sorry("Resolution from half-maps already worse than target")
          else:
            fsc_results = wmc1.fsc(other=wmc2)
            start_FSC = fsc_results.fsc
            start_d_inv = fsc_results.d_inv
        else:
          mc2 = mc1.deep_copy()
          wmc2 = mc2.deep_copy()

        target_params = fsc_params_from_d_min(target_d_min)
        r_target = target_params.signalratio
        dB_target = target_params.deltaB
        nref_total = mc1.size()
        nref_bin1 = 50. # target number of reflections in first bin
        nbins = max(6, round(nref_total / nref_bin1)**(2./3.))
        wmc1.setup_binner_d_star_sq_bin_size(max_bins=nbins)
        wmc2.use_binner_of(wmc1)
        n_too_poor = 0 # Keep track of how many terms are worse than target
        for i_bin in wmc1.binner().range_used():
          sel = wmc1.binner().selection(i_bin)
          wmc1sel = wmc1.select(sel)
          wmc2sel = wmc2.select(sel)
          nrefsel = wmc1sel.data().size()
          mean_ssqr = flex.mean(wmc1sel.d_star_sq().data())
          # d_bin = math.sqrt(1./mean_ssqr)
          if halfmaps:
            sigmaT = flex.mean(flex.pow2(wmc1sel.amplitudes().data()) +
                              flex.pow2(wmc2sel.amplitudes().data())) / 2
            sigmaE_0 = flex.mean(flex.pow2(flex.abs(wmc1sel.data() - wmc2sel.data())) / 2)
            sigmaS = max(0., sigmaT - sigmaE_0) # Negative correlation -> zero signal
          else: # Approximate input map as perfect to target resolution
            sigmaS = flex.mean(flex.pow2(wmc1sel.amplitudes().data()))
            sigmaE_0 = 0.
          sigmaE_factor = math.exp(dB_target * mean_ssqr / 4.) / r_target
          sigmaE_factor = min(1000., sigmaE_factor) # Avoid overkill on noise
          sigmaE_target = sigmaS * sigmaE_factor
          delta_sigmaE = max(0., sigmaE_target - sigmaE_0) # Can't subtract noise
          # This delta_sigmaE would account for noise in masked region alone
          delta_sigmaE /= mask_fraction # Correct for spreading over full box
          if delta_sigmaE > 0.:
            target_sigma = math.sqrt(delta_sigmaE/2.) # For real and imaginary parts
            random_complex = flex.complex_double(
                                np.random.normal(scale=target_sigma,size=nrefsel) +
                            1j*np.random.normal(scale=target_sigma,size=nrefsel))

            mc1_new.data().set_selected(sel, mc1.data().select(sel) + random_complex )
            random_complex = flex.complex_double(
                                np.random.normal(scale=target_sigma,size=nrefsel) +
                            1j*np.random.normal(scale=target_sigma,size=nrefsel))
            mc2_new.data().set_selected(sel, mc2.data().select(sel) + random_complex )
          else:
            mc1_new.data().set_selected(sel, mc1.data().select(sel))
            mc2_new.data().set_selected(sel, mc2.data().select(sel))
            n_too_poor += nrefsel

        frac_too_poor = n_too_poor / mc1_new.data().size()
        if frac_too_poor > 0.25:
          raise Sorry("Initial resolution seems to be worse than target")
        mmm.add_map_from_fourier_coefficients(mc1_new, map_id='map_manager_1')
        mmm.add_map_from_fourier_coefficients(mc2_new, map_id='map_manager_2')
        working_mmm = apply_spherical_mask(mmm, target_d_min)
        wmc1 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_1')
        wmc2 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=1000., map_id='map_manager_2')
        d_min_results = d_min_from_fsc_analytical(wmc1, wmc2, guess_d_min=target_d_min)
        d_min_from_fsc = d_min_results.d_min
        end_d_inv = d_min_results.d_inv
        end_FSC = d_min_results.fsc_bins

        # Make output half-maps
        filenames = []
        output_dir = params.output_files.output_dir

        if output_dir and output_dir != ".":
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            output_root = os.path.join(output_dir, file_root)
        else:
            output_root = file_root
        # Use data_manager for writing
        dm.set_overwrite(True)
        mc1_new = mc1_new.select(mc1_new.d_spacings().data() >= target_d_min)
        mmm.add_map_from_fourier_coefficients(mc1_new, map_id='map1')
        mm1 = mmm.get_map_manager_by_id(map_id='map1')
        output_path1 = output_root + '_half_map_1.map'
        dm.write_real_map_file(mm1, filename=output_path1)
        filenames.append(output_path1)
        mc2_new = mc2_new.select(mc2_new.d_spacings().data() >= target_d_min)
        mmm.add_map_from_fourier_coefficients(mc2_new, map_id='map2')
        mm2 = mmm.get_map_manager_by_id(map_id='map2')
        output_path2 = output_root + '_half_map_2.map'
        dm.write_real_map_file(mm2, filename=output_path2)
        filenames.append(output_path2)

        # Average and write out reduced resolution half-maps
        mm_mean_data = (mm1.map_data() + mm2.map_data()) / 2
        mm_mean = mm1.customized_copy(map_data=mm_mean_data)
        output_pathf = output_root + '_map.map'
        dm.write_real_map_file(mm_mean, filename=output_pathf)
        filenames.append(output_pathf)

        if not mute:
          print("Simulated data d_min_from_fsc: ", d_min_from_fsc, file=log)
          print("Reduced resolution map files written to:\n",
                "    " + output_path1 + "\n",
                "    " + output_path2 + "\n",
                "    " + output_pathf)

        if plot:
          import matplotlib.pyplot as plt
          analytical_FSC = []
          for id in range(len(end_d_inv)):
            d_inv = end_d_inv[id]
            analytical_FSC.append(FSCfunction(r_target, dB_target, d_inv))
          fig, ax = plt.subplots(1,1)
          ax.set_xlabel('1/d')
          ax.set_ylabel('FSC')
          if halfmaps:
            ax.plot(start_d_inv, start_FSC, label="Starting FSC")
          ax.plot(end_d_inv, end_FSC, label="Final FSC")
          ax.plot(end_d_inv, analytical_FSC, label="Target FSC")
          ax.axvline(x=1./target_d_min, color="grey", linestyle=":",
                     label="Filter resolution")
          ax.legend()
          output_path = os.path.join(output_root + "_FSC.png")
          plt.savefig(output_path)
          filenames.append(output_path)
          if not mute:
            print("FSC curves written to: " + output_path)

        self.d_min_from_fsc = d_min_from_fsc
        self.filenames = filenames

    def get_results(self):
        return group_args(d_min_from_fsc=self.d_min_from_fsc,
                          filenames = self.filenames)


# -----------------------------------------------------------------------------
# Standard Entry Point
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    run_program(program_class=Program)
