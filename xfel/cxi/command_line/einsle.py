from __future__ import division
from __future__ import print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME LM14.einsle
import sys,os
import iotbx.pdb
import string
from scitbx.array_family import flex
import mmtbx.command_line.fmodel
import mmtbx.f_model

"""Fit the anomalous scattering parameters f' and f", given the pdb model (for phases) and the
   intensities with Friedel mates separated, as in Oliver Einsle, Susana L. A. Andrade, Holger Dobbek,
   Jacques Meyer, and Douglas C. Rees (2007). Assignment of Individual Metal Redox States in a
   Metalloprotein by Crystallographic Refinement at Multiple X-ray Wavelengths. Journal of the
   American Chemical Society 129, 2210-2211."""

class Model:
  def __init__(self,file_name, d_min, algorithm = "direct", use_solvent=False, plot=True):
    """Can we compare Fmodel calculated from Phenix GUI, vs. from a script?
       Answer: they are in perfect agreement.
       if False: script calculated f_calc without solvent
       if True:  script calculated f_model with solvent. OK"""
    self.d_min = d_min
    self.pdb_inp = iotbx.pdb.input(file_name)
    self.xray_structure = self.pdb_inp.xray_structure_simple()
    self.xray_structure.show_summary()
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    params2.output.type = "complex"
    params2.high_resolution = d_min
    if use_solvent :
      params2.fmodel.k_sol = 0.35
      params2.fmodel.b_sol = 46.
    self.f_model_complex = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = params2).f_model

    self.f_model_real = abs(self.f_model_complex)
    self.f_model_real.set_observation_type_xray_amplitude()
    self.f_model_real.show_summary(prefix="FMODEL ")

  def wavelength_independent_phases(self,group_sulfurs):
    complex = self.f_model_complex
    self.phases=complex.phases(deg=False)

    #now loop through the heavy atom scatterers
    pdb_hierarchy = self.pdb_inp.construct_hierarchy()
    asc = pdb_hierarchy.atom_selection_cache()
    selection = asc.selection("(element Zn or element Ca or element S or element Fe or element Yb)")
    #selection = asc.selection("element Zn or element Fe")
    #selection = asc.selection("chain A and resseq 201")
    print("%d atoms selected out of total %d"%(
    selection.count(True), selection.size()))
    self.N_anom_scatterers = selection.count(True)
    xrs_ions    = self.xray_structure.select(selection)
    xrs_ions.show_scatterers()
    self.scatterer_model_idx = self.get_scatterer_model_idx(xrs_ions,group_sulfurs)
    self.scatterer_idx = selection.iselection(True)
    self.f_calc = self.xray_structure.structure_factors(d_min=self.d_min, algorithm=algorithm).f_calc()
    self.metals = xrs_ions

  def get_scatterer_model_idx(self,ions,group_sulfurs):
    if not group_sulfurs:
      return flex.size_t(range(len(ions.scatterers())))
    # if it is desired to group all the sulfurs together
    sulfur_idx = None
    parameter_no = 0
    indirection = flex.size_t()
    for j,sc in enumerate(ions.scatterers()):
      if sc.scattering_type=="S":
        if sulfur_idx is None:
          sulfur_idx = parameter_no
          parameter_no += 1
        indirection.append(sulfur_idx)
      else:
        indirection.append(parameter_no)
        parameter_no+=1
    return indirection

  def parts(self,fpp):
    dano_summation = None
    Ddanocalc_Dp = []
    for j,idx in enumerate(self.scatterer_idx):
      bool_array = self.xray_structure.by_index_selection([idx])
      xrs_atom = self.xray_structure.select(bool_array)
      #xrs_atom.show_scatterers()
      f_calc_atom = self.f_calc.structure_factors_from_scatterers(
        xray_structure=xrs_atom, algorithm=algorithm).f_calc()
      adata = f_calc_atom.data().parts()[0]
      bdata = f_calc_atom.data().parts()[1]

      product_factor = bdata * flex.cos(self.phases.data()) - adata * flex.sin(self.phases.data())
      # correspnds to -b cos(alpha) + a sin(alpha)
      #print list(xrs_atom.scattering_types())
      #xrs_atom.scattering_type_registry().show_summary()
      #xrs_atom.scattering_type_registry().show()
      f_zero = xrs_atom.scattering_type_registry().sum_of_scattering_factors_at_diffraction_angle_0()
      #print f_zero
      Ddanocalc_Dp.append((-2./f_zero)*product_factor)
      term = (-2.*fpp[j]/f_zero)*product_factor
      if dano_summation is None:
        dano_summation = term
      else:
        dano_summation += term
    return dano_summation,Ddanocalc_Dp


    """
Next things to do:
DONE Evaluate the f" summation, looping over all the heavy atom scatterers
DONE Calculate functional and gradient
DONE minimize and inspect results for thermolysin
DONE do the same for the ferredoxin
do the same for the high-energy remote of ferredoxin
DONE Refine the LD91 lysozyme model
DONE try f" refinement on lysozyme

debugs:
DONE Report variance-weighted C.C.
DONE Wrap all Sulfurs into a single parameter
Fit the f' as well
---->Fit Obs to the Fcalc with a B-factor or bin scaling
DONE Try to sort out f" for a nested subset of Yb-lyso datasets.
DONE Compare my f" refined values with tabular values (in fact, look them up with a script)
-->make sure Mona's program outputs the wavelength
---->Sort out the geometry in AxFd monomer A
-->Look at the high-energy remote of ferredoxin
Phil-out a few cases so we can get the pdb files out of the code
"""

  def scale(self,other):
    from cctbx import miller
    matches = miller.match_indices(self.f_model_real.indices(),other.indices())
    sel0 = flex.size_t([p[0] for p in matches.pairs()])
    sel1 = flex.size_t([p[1] for p in matches.pairs()])

    val0 = self.f_model_real.data().select(sel0)
    val1 = other.data().select(sel1)
    plot=False
    if plot:
      from matplotlib import pyplot as plt
      plt.plot([-1,4],[-1,4],"g-")
      plt.plot(flex.log10(val0),flex.log10(val1),"r.")
      plt.show()

    from xfel.cxi.cxi_cc import correlation
    slope,offset,corr,N = correlation(
      self = self.f_model_real.select(sel0),
      other = other.select(sel1))
    print(slope,offset,corr,N)
    if plot:
      from matplotlib import pyplot as plt
      plt.plot([-1,4],[-1,4],"g-")
      plt.plot(flex.log10(val0),flex.log10(slope * val1),"r,")
      plt.show()
    return slope

  def scaling_metrics(self,other):
    # Read reflections
    # some requirements. 1) Fobs scaled to Fcalc, not the other way around.
    # 2) ability to make a plot of the two scaled sets
    # 3) set the number of bins
    # 4) understand and print out the per-bin scaling factor
    # 5) print an overall stats line at the end of the table
    # 6) choose one or the other binnings
    """1) scaling and analysis are separate functions"""

    #f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
    f_obs = other
    #r_free_flags = r_free_flags.array(data=r_free_flags.data()==1)
    # Read model

    # Get Fmodel
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      #r_free_flags   = r_free_flags,
      xray_structure = self.xray_structure)
    # Do anisotropic overall scaling, bulk-solvent modeling, outlier rejection
    #fmodel.update_all_scales()
    print("r_work, r_free: %6.4f, %6.4f"%(fmodel.r_work(), fmodel.r_free()))
    # Print statistics in resolution bins
    f_model = fmodel.f_model_scaled_with_k1()
    bin_selections = fmodel.f_obs().log_binning()
    dsd = fmodel.f_obs().d_spacings().data()
    print("Bin# Resolution    Nref Cmpl  Rw     CC")
    fmt="%2d: %6.3f-%-6.3f %5d %5.3f %6.4f %6.4f"
    for i_bin, sel in enumerate(bin_selections):
      d           = dsd.select(sel)
      d_min       = flex.min(d)
      d_max       = flex.max(d)
      fmodel_sel  = fmodel.select(sel)
      n           = d.size()
      f_obs_sel   = fmodel.f_obs().select(sel)
      f_model_sel = abs(f_model.select(sel)).data()
      cmpl        = f_obs_sel.completeness(d_max=d_max)
      r_work      = fmodel_sel.r_work()
      cc          = flex.linear_correlation(x=f_obs_sel.data(),
                    y=f_model_sel).coefficient()
      print(fmt%(i_bin, d_max, d_min, n, cmpl, r_work, cc))
    # Alternative binning
    print()
    print("Bin# Resolution    Nref Cmpl  Rw     CC")
    fmodel.f_obs().setup_binner(reflections_per_bin = 2500)
    f_model.use_binning_of(fmodel.f_obs())
    for i_bin in fmodel.f_obs().binner().range_used():
      sel = fmodel.f_obs().binner().selection(i_bin)
      d           = dsd.select(sel)
      d_min       = flex.min(d)
      d_max       = flex.max(d)
      fmodel_sel  = fmodel.select(sel)
      n           = d.size()
      f_obs_sel   = fmodel.f_obs().select(sel)
      f_model_sel = abs(f_model.select(sel)).data()
      cmpl        = f_obs_sel.completeness(d_max=d_max)
      r_work      = fmodel_sel.r_work()
      cc          = flex.linear_correlation(x=f_obs_sel.data(),
                    y=f_model_sel).coefficient()
      print(fmt%(i_bin, d_max, d_min, n, cmpl, r_work, cc))


def get_obs(file_name,tag):
  """Can we scale one amplitude array to another?"""
  Possible=["i(+)","iobs(+)"]
  if tag is not None:  Possible.append(tag.lower())
  from iotbx import mtz
  data_SR = mtz.object(file_name)
  for array in data_SR.as_miller_arrays():
    this_label = array.info().label_string().lower()
    array.show_summary(prefix="OBS ")
    if True in [this_label.find(tag)>=0 for tag in Possible]: break
  assert True in [this_label.find(tag)>=0 for tag in Possible], \
         "Cannot find i(+); use phenix.mtz.dump and give iobs_tag in phil string"
  wavelength = get_wavelength(data_SR,Possible)
  f_ampl = array.as_amplitude_array()
  merged = array.average_bijvoet_mates()
  f_ampl_merged = merged.as_amplitude_array()

  return f_ampl,f_ampl_merged,wavelength

def get_wavelength(self, Possible):#mtz_obj
    # code is adapted from iotbx/mtx/__init__.py show_summary()
    wavelength = None
    for i_crystal,crystal in enumerate(self.crystals()):
      for i_dataset,dataset in enumerate(crystal.datasets()):
        if (dataset.n_columns() > 0):
          fields_list = [[
            "label", "#valid", "%valid", "min", "max", "type", ""]]
          max_field_lengths = [len(field) for field in fields_list[0]]
          max_field_lengths[-2] = 0
          for i_column,column in enumerate(dataset.columns()):
            fields = column.format_fields_for_mtz_dump(
              n_refl=self.n_reflections())
            fields_list.append(fields)
            for i,field in enumerate(fields):
              max_field_lengths[i] = max(max_field_lengths[i], len(field))
          format = "    %%-%ds %%%ds %%%ds %%%ds %%%ds %%%ds %%s" % tuple(
            max_field_lengths[:6])
          for fields in fields_list:
            if fields[0].lower().strip() in Possible:
              wavelength = dataset.wavelength()
    return wavelength


import scitbx.lbfgs
class FPP_optimizer:
  def __init__(self, model,diffs,params):
    self.params = params
    self.model = model
    if params.group_sulfurs:
      self.n = 1 + flex.max(self.model.scatterer_model_idx)
    else:
      self.n = self.model.N_anom_scatterers
    self.x = flex.double(self.n,0.)

    from cctbx import miller
    matches = miller.match_indices(self.model.f_model_real.indices(),diffs.indices())
    self.sel0 = flex.size_t([p[0] for p in matches.pairs()])
    self.sel1 = flex.size_t([p[1] for p in matches.pairs()])

    self.diffs = diffs.select(self.sel1)

    print("SELECTED %d diffs out of %d"%(len(self.diffs.data()), len(diffs.data())))

    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-4,
        max_iterations=20))

  def compute_functional_and_gradients(self,plot=False):
    dano_summation,Ddanocalc_Dp = self.model.parts(self.get_fpp())
    if plot: self.plot(dano_summation)

    #calculate the functional
    residual = self.diffs.data() - dano_summation.select(self.sel0)
    if self.params.use_weights:  residual /= self.diffs.sigmas()
    F = 0.5 * flex.sum(residual * residual)
    print(("LBFGS stp",F))

    g = flex.double(self.n, 0.)
    for j,item in enumerate(Ddanocalc_Dp):
      deriv = item.select(self.sel0)
      vector=-deriv*residual
      if self.params.use_weights:  vector /= self.diffs.sigmas()

      g[self.model.scatterer_model_idx[j]] += flex.sum(vector)
    return F, g

  def correlation(self):
    dano_summation,Ddanocalc_Dp = self.model.parts(self.get_fpp())

    #calculate the correlation
    self.diffs.data(),dano_summation.select(self.sel0)
    if self.params.use_weights:
      wt = 1./(self.diffs.sigmas()*self.diffs.sigmas())
    else: wt = flex.double(len(self.sel0), 1.)

    from scitbx.math.tests.tst_weighted_correlation import weighted_correlation
    return weighted_correlation(wt, self.diffs.data(), dano_summation.select(self.sel0))

  def get_fpp(self):
    if self.params.group_sulfurs:
      return self.x.select(self.model.scatterer_model_idx)
    else:
      return self.x

  def plot(self,dano_summation):
    from matplotlib import pyplot as plt

    if self.params.use_weights:
      wt = 1./(self.diffs.sigmas()*self.diffs.sigmas())
      order = flex.sort_permutation(wt)
      wt = wt.select(order)
      df = self.diffs.data().select(order)
      dano = dano_summation.select(self.sel0).select(order)
      from matplotlib.colors import Normalize
      dnorm = Normalize()
      dnorm.autoscale(wt.as_numpy_array())
      CMAP = plt.get_cmap("rainbow")
      for ij in range(len(self.diffs.data())):
        #blue represents zero weight:  red, large weight
        plt.plot([df[ij]],[dano[ij]],color=CMAP(dnorm(wt[ij])),marker=".", markersize=4)

    else:
      plt.plot(self.diffs.data(),dano_summation.select(self.sel0),"r,")
    plt.axes().set_aspect("equal")
    plt.axes().set_xlabel("Observed Dano")
    plt.axes().set_ylabel("Model Dano")
    plt.show()


def show_scatterers(self,fpp,wave):
    print ("""Label            f"      Coordinates        Occupancy """
                 "Uiso, Ustar as Uiso")
    scatterers = self.scatterers()
    types_used = {}
    for j,sc in enumerate(scatterers):
      sc.fdp = fpp[j]
      show_scatterer(sc, unit_cell=self.unit_cell())
      types_used[sc.scattering_type]=None
    print("Tabular values for %7.1f eV (%8.5f Angstrom):"%(12398/wave,wave))
    from cctbx.eltbx import sasaki, henke
    for tag in types_used.keys():
      fpp_expected_sasaki = sasaki.table(tag).at_angstrom(
          wave).fdp()
      fpp_expected_henke = henke.table(tag).at_angstrom(
          wave).fdp()
      print("           %-4s"%tag, end=' ')
      print("%6.3f" % (max(fpp_expected_sasaki,fpp_expected_henke)))

def show_scatterer(self, f=None, unit_cell=None):

  if (f is None): f = sys.stdout
  from cctbx import adptbx
  print("%-4s" % self.label[5:15], end=' ', file=f)
  print("%-4s" % self.scattering_type, end=' ', file=f)
  #print >> f, "%3d" % self.multiplicity(),
  print("%6.3f" % (self.fdp), end=' ', file=f)
  print("(%7.4f %7.4f %7.4f)" % self.site, end=' ', file=f)
  print("%4.2f" % self.occupancy, end=' ', file=f)
  if self.flags.use_u_iso():
    print("%6.4f" % self.u_iso, end=' ', file=f)
  else:
    print('[ - ]', end=' ', file=f)
  if self.flags.use_u_aniso():
    assert unit_cell is not None
    u_cart = adptbx.u_star_as_u_cart(unit_cell, self.u_star)
    print("%6.4f" % adptbx.u_cart_as_u_iso(u_cart), file=f)
    print("     u_cart =", ("%6.3f " * 5 + "%6.3f") % u_cart, end=' ', file=f)
  else:
    print('[ - ]', end=' ', file=f)
  if False and (self.fp != 0 or self.fdp != 0):
    print("\n     fp,fdp = %6.4f,%6.4f" % (
      self.fp,
      self.fdp), end=' ', file=f)
  print(file=f)

from libtbx.phil import parse
phil_scope = parse("""
  pdb = None
    .type = str
  mtz = None
    .type = str
  d_min = 2.0
    .type = float
  use_weights = True
    .type = bool
    .help = Use variance weighting for LSQ fitting and plotting
  make_plot = True
    .type = bool
    .help = Plot everything at the end
  group_sulfurs = True
    .type = bool
    .help = Treat all sulfur atoms as having the same f' and f"
  iobs_tag = None
    .type = str
    .help = lower-case string identifying the anomalous iobs(+)
    .help = need not give i(+) or iobs(+) as these are already hard coded
  wavelength_overrride = None
    .type = float
  """, process_includes = True)

def get_params(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      try:
        user_phil.append(parse(file_name=arg))
      except Exception as e:
        print(str(e))
        raise Sorry("Couldn't parse phil file %s"%arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        print(str(e))
        raise Sorry("Couldn't parse argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params


if __name__ == "__main__":
  params = get_params(sys.argv[1:])
  #params.mtz = "/Users/nksauter/xtalwork/LM14/Yb-lysozyme/merge118/Yblyso118_post_anom_4etc_s0_mark0.mtz"
  params.wavelength_overrride=1.3853 # for Yb-lysozyme
  obs_ampl,obs_ampl_merged,wave = get_obs(params.mtz, params.iobs_tag)
  if wave==None or wave==1.:  wave = params.wavelength_overrride
  if params.d_min == None:
    params.d_min = obs_ampl_merged.d_min()
  print("DMIN = %7.2f"%params.d_min)

  #params.pdb = "/Users/nksauter/xtalwork/LM14/Yb-lysozyme/Refine_4/LM14_1colorYblyso_refine_4.pdb"
  algorithm = "direct"
  algorithm = "fft"
  M = Model(params.pdb,params.d_min)
  slope = M.scale(other=obs_ampl_merged)

  #Note the observations are scaled to the Fmodel, not the other way around
  scaled_obs_ampl = obs_ampl.customized_copy(data = slope*obs_ampl.data(), sigmas=slope*obs_ampl.sigmas())
  scaled_obs_ampl_merged = obs_ampl_merged.customized_copy(
    data = slope*obs_ampl_merged.data(), sigmas=slope*obs_ampl_merged.sigmas())

  M.scaling_metrics(other = scaled_obs_ampl)

  anomalous_diffs = scaled_obs_ampl.anomalous_differences()

  for x in range(5):
    print(anomalous_diffs.indices()[x], anomalous_diffs.data()[x],anomalous_diffs.sigmas()[x])

  M.wavelength_independent_phases(params.group_sulfurs)
  print(M.N_anom_scatterers ,"anomalous_scatterers")
  M.parts(flex.double(M.N_anom_scatterers))
  FPPO = FPP_optimizer(M, anomalous_diffs, params)

  show_scatterers(M.metals, FPPO.get_fpp(), wave)
  print("C.C. = %7.2f%%"%(100. * FPPO.correlation()))
  if params.make_plot:
    FPPO.compute_functional_and_gradients(plot=True)
