
from __future__ import division
from mmtbx.validation import residue, atom, validation
from libtbx import slots_getstate_setstate
from libtbx.utils import null_out
import sys

__real_space_attr__ = [
  "b_iso",
  "fofc",
  "two_fofc",
  "fmodel",
]

class residue_real_space (residue) :
  # CC is 'score' attribute
  __slots__ = residue.__slots__ + __real_space_attr__

  @property
  def cc (self) :
    return self.score

  @staticmethod
  def header () :
    return "%-20s  %6s  %4s  %6s  %6s  %5s" % ("atom", "b_iso", "occ",
      "2Fo-Fc", "Fmodel", "CC")

  def as_string (self) :
    return "%-20s  %6.2f  %4.2f  %6.2f  %6.2f  %5.3f" % (self.id_str(),
      self.b_iso, self.occupancy, self.two_fofc, self.fmodel, self.score)

class data_statistics (slots_getstate_setstate) :
  __slots__ = [
    "d_max",
    "d_min",
    "info",
    "n_refl",
    "r_work",
    "r_free",
    "twin_law",
    "wilson_b",
    "r_work_outer",
    "r_free_outer",
    "d_max_outer",
    "d_min_outer",
    "completeness_outer",
    "n_refl_outer",
  ]
  def __init__ (self, fmodel) :
    f_obs = fmodel.f_obs().deep_copy()
    f_obs.setup_binner(n_bins=10)
    self.d_max = f_obs.d_max_min()[0]
    self.d_min = f_obs.d_min()
    self.info = fmodel.info()
    self.n_refl = f_obs.indices().size()
    self.r_free = self.info.r_free
    self.r_work = self.info.r_work
    self.twin_law = fmodel.twin_law
    self.wilson_b = fmodel.wilson_b()
    # outer shell
    d_max_min_outer = f_obs.binner().bin_d_range(10)
    self.d_max_outer = d_max_min_outer[0]
    self.d_min_outer = d_max_min_outer[1]
    self.r_free_outer = fmodel.r_free(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    self.r_work_outer = fmodel.r_work(d_max=self.d_max_outer,
      d_min=self.d_min_outer)
    self.completeness_outer = f_obs.completeness(d_max=self.d_max_outer)
    self.n_refl_outer = f_obs.resolution_filter(d_max=self.d_max_outer,
      d_min=self.d_min_outer).indices().size()

  def show_summary (self, out=sys.stdout, prefix="") :
    print >> out, "%sHigh resolution       = %7.3f" % (prefix, self.d_min)
    print >> out, "%sR-work                = %8.4f" % (prefix, self.r_work)
    print >> out, "%sR-free                = %8.4f" % (prefix, self.r_free)

  def show (self, out=sys.stdout, prefix="") :
    print >> out, "%sResolution range      = %7.3f - %.3f (%.3f - %.3f)" % \
      (prefix, self.d_max, self.d_min, self.d_max_outer, self.d_min_outer)
    print >> out, "%sNumber of reflections = %8d (%d)" % (prefix, self.n_refl,
      self.n_refl_outer)
    print >> out, "%sCompleteness          = %6.2f%% (%.2f%%)" % (prefix,
      self.info.completeness_in_range*100, self.completeness_outer*100)
    print >> out, "%sR-work                = %8.4f (%.4f)" % (prefix,
      self.r_work, self.r_work_outer)
    print >> out, "%sR-free                = %8.4f (%.4f)" % (prefix,
      self.r_free, self.r_free_outer)

class real_space (validation) :
  """
  Real-space correlation calculation for residues with at least two atoms.
  """
  def get_result_class (self) : return residue_real_space

  def __init__ (self, fmodel, pdb_hierarchy, cc_min=0.8) :
    validation.__init__(self)
    from mmtbx import real_space_correlation
    try :
      rsc_params = real_space_correlation.master_params().extract()
      rsc_params.detail="residue"
      rsc_params.map_1.fill_missing_reflections = False
      rsc_params.map_2.fill_missing_reflections = False
      rsc = real_space_correlation.simple(
        fmodel=fmodel,
        pdb_hierarchy=pdb_hierarchy,
        params=rsc_params,
        log=null_out())
    except Exception, e :
      print >> out, "Error: %s" % str(e)
    else :
      rsc_by_res = []
      for i, result_ in enumerate(rsc) :
        if (result_.n_atoms == 1) or (result_.residue.resname == "HOH") :
          continue
        result = residue_real_space(
          chain_id=result_.chain_id,
          resname=result_.residue.resname,
          resseq=result_.residue.resseq,
          icode=result_.residue.icode,
          altloc="",
          score=result_.cc,
          b_iso=result_.b,
          occupancy=result_.occupancy,
          fmodel=result_.map_value_1,
          two_fofc=result_.map_value_2,
          outlier=result_.cc < cc_min,
          xyz=result_.residue.atoms().extract_xyz().mean())
        if result.is_outlier() :
          self.n_outliers += 1
        # XXX unlike other validation metrics, we always save the results for
        # the real-space correlation, since these are used as the basis for
        # the multi-criterion plot in Phenix.  The show() method will only
        # print outliers, however.
        self.results.append(result)

  def show_summary (self, out=sys.stdout, prefix="") :
    print >> out, prefix + "%d residues with CC(Fc,2mFo-DFc) < 0.8" % \
      self.n_outliers

  def show (self, out=sys.stdout, prefix="  ", verbose=True) :
    if (self.n_outliers > 0) :
      print >> out, prefix + self.get_result_class().header()
      for result in self.results :
        if result.is_outlier() :
          print >> out, prefix + str(result)
    self.show_summary(out=out, prefix=prefix)
