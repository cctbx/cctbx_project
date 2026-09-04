from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.f_model
import mmtbx.bulk_solvent.scaler
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import iotbx.phil

llgi_bss_params = iotbx.phil.parse("""\
  enabled = False
    .type = bool
    .short_caption = Fit bulk-solvent/scaling against Feff instead of F-obs
    .help = "If True (and target=llgi), redo the whole least-squares "\
            "bulk-solvent-and-scaling fit (overall scale, anisotropic "\
            "scaling, and, unless llgi_data.e_scale_bulk_solvent.enabled "\
            "is also True, k_sol/B_sol) against llgi_data.feff instead "\
            "of fmodel.f_obs(). Runs as a distinct task immediately "\
            "after the ordinary (F-obs) bss task, overwriting its "\
            "k_isotropic/k_anisotropic/k_mask result -- see "\
            "run_llgi_bss. fmodel.f_obs() itself is never modified, so "\
            "R-work/R-free and map coefficients (still F-obs-based, see "\
            "mmtbx.f_model.manager._r_factor/electron_density_map) are "\
            "UNCHANGED by this option; only the scale/bulk-solvent "\
            "parameters feeding sigmaA/ScatFrac estimation and the LLGI "\
            "coordinate-refinement target itself are affected. This is "\
            "still an ordinary least-squares fit (same machinery as the "\
            "F-obs bss step, mmtbx.bulk_solvent.scaler), not an LLGI-"\
            "likelihood fit -- see llgi_e_bulk_solvent.py's E-scale "\
            "machinery for that (a separate, independent option, whose "\
            "own k_sol/B_sol result -- if enabled -- takes precedence "\
            "over this task's, since it runs afterward in the "\
            "macrocycle task order). Off by default: new, untested "\
            "machinery layered on top of the already-DEVELOPER_ONLY llgi "\
            "target."
  fast = True
    .type = bool
    .expert_level = 3
    .help = "Use the fast (analytical) scaler, matching bss's own "\
            "default (update_all_scales(fast=True)), rather than the "\
            "slower minimization-based one."
""")

def run_llgi_bss(fmodel, params=None, log=None):
  """ Redo the least-squares bulk-solvent-and-scaling fit (overall scale,
  anisotropic scaling, and bulk solvent -- k_sol/B_sol, via a single
  exponential shell, matching bss's own fast/analytical default) against
  llgi_data.feff instead of fmodel.f_obs(), and apply the result
  (k_isotropic, k_anisotropic, k_mask) back onto fmodel in place.

  Meant to run as a distinct task immediately after the ordinary F-obs
  bss task (see phenix.refinement.macro_cycle.py's updatellgibss, which
  mirrors updatellgiebulksolvent's own placement/gating pattern) --
  bss has already produced a scaled fmodel (overall scale + anisotropy +
  bulk solvent, all fit against F-obs), which supplies this function's
  starting point (k_sol/b_sol carried over from bss's own core, matching
  llgi_e_bulk_solvent.run_inner_loop's own "fmodel already scaled"
  precondition) and, via fmodel.f_calc()/f_masks(), the same structure-
  factor arrays bss itself just used. This function then REPLACES
  k_isotropic/k_anisotropic/k_mask with the LS fit against Feff.

  Unlike llgi_e_bulk_solvent.run_inner_loop (which fits sigmaA(d) and
  bulk solvent by alternating LLGI-likelihood optimisation, and
  deliberately leaves overall scale/anisotropy exactly as bss set them),
  this is an ordinary least-squares fit -- the same scaler machinery bss
  itself uses (mmtbx.bulk_solvent.scaler.run_simple /
  mmtbx.bulk_solvent.bulk_solvent_and_scaling.bulk_solvent_and_scales),
  just pointed at a different amplitude array. Deliberately simple: no
  outlier rejection here (nacelle's own upstream INFO/EINFO-based
  reflection selection -- see phenix.refinement.llgi_data.get_llgi_data
  -- already plays that role for Feff; re-running F-obs-style amplitude/
  sigma outlier rejection against Feff, which has no meaningful sigmas,
  would not be well-founded -- see doc/llgi_target_design.md's Feff bss
  note), and no riding-hydrogen scattering re-optimisation (left exactly
  as bss's own F-obs fit set it -- H scattering is a small, F-obs-scale-
  insensitive correction, not worth refitting here).

  fmodel: an mmtbx.f_model.manager, already scaled by the ordinary F-obs
  bss task (see above) and carrying llgi_data (feff/dobs/teps/resn,
  already matched to fmodel.f_obs()'s current index set -- see
  phenix.refinement.llgi_data.get_llgi_data). Raises Sorry if llgi_data
  is not attached.

  params: extracted llgi_bss_params phil, or None for defaults.

  Returns a group_args with .k_sol, .b_sol, .b_cart (the fitted bulk-
  solvent/anisotropy parameters, as reported by mmtbx.bulk_solvent.
  scaler.run for logging/diagnostics -- same shape as f_model_all_scales.
  run's own return value) and .r_all_feff (the LS scaler's own R-factor
  against Feff, NOT comparable to fmodel.r_work()/r_free(), which remain
  F-obs-based throughout -- see this function's docstring above).
  """
  from libtbx.utils import Sorry
  if(params is None):
    params = llgi_bss_params.extract()
  llgi_data = fmodel.llgi_data()
  if(llgi_data is None):
    raise Sorry(
      "run_llgi_bss requires llgi_data (FEFF) to already be attached "
      "(see mmtbx.f_model.manager.set_llgi_data).")
  feff = llgi_data.feff
  if(not feff.indices().all_eq(fmodel.f_obs().indices())):
    raise Sorry(
      "run_llgi_bss: llgi_data.feff's index set does not match "
      "fmodel.f_obs()'s current index set (has f_obs shrunk, e.g. via "
      "outlier removal, since llgi_data was last attached/refreshed?).")
  fmodel_kbu = mmtbx.f_model.manager_kbu(
    f_obs   = feff,
    f_calc  = fmodel.f_calc(),
    f_masks = fmodel.f_masks(),
    ss      = fmodel.ss)
  bss_params = bss.master_params.extract()
  if(params.fast):
    result = mmtbx.bulk_solvent.scaler.run_simple(
      fmodel_kbu     = fmodel_kbu,
      r_free_flags   = fmodel.r_free_flags(),
      bulk_solvent   = bss_params.bulk_solvent,
      aniso_scale    = bss_params.anisotropic_scaling,
      bin_selections = fmodel.bin_selections)
    k_sol, b_sol = None, None
  else:
    result = bss.bulk_solvent_and_scales(
      fmodel_kbu = fmodel_kbu,
      params     = bss_params)
    k_sol = result.k_sols()
    b_sol = result.b_sols()
  r_all_feff = result.r_all()
  fmodel.update_core(
    k_mask        = result.k_masks(),
    k_anisotropic = result.k_anisotropic(),
    k_isotropic   = result.k_isotropic())
  if(log is not None):
    print(
      "  LLGI bss (Feff-based): R-all(Feff)=%.4f  r_work(F-obs)=%.4f  "
      "r_free(F-obs)=%.4f" % (r_all_feff, fmodel.r_work(), fmodel.r_free()),
      file=log)
  from libtbx import group_args
  return group_args(k_sol=k_sol, b_sol=b_sol, r_all_feff=r_all_feff)
