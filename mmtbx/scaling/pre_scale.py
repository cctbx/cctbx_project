from cctbx import miller
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews
from mmtbx.scaling import data_statistics


class pre_scaler(object):
  def __init__(self,
               miller_array,
               pre_scaling_protocol,
               basic_info,
               out=None):
    ## Make deep copy of the miller array of interest
    self.x1 = miller_array.deep_copy()
    self.options=pre_scaling_protocol
    self.basic_info= basic_info

    ## Determine unit_cell contents
    print >> out
    print >> out, "Matthews analyses"
    print >> out, "-----------------"
    print >> out
    print >> out, "Inspired by: Kantardjieff and Rupp. Prot. Sci. 12(9): 1865-1871 (2003)."
    matthews_analyses = matthews.matthews_rupp(
      miller_array = self.x1,
      n_residues = self.basic_info.n_residues,
      n_bases = self.basic_info.n_bases,
      out=out, verbose=1)
    n_residues=matthews_analyses[0]
    n_bases=matthews_analyses[1]
    n_copies_solc=matthews_analyses[2]

    if (self.basic_info.n_residues==None):
      self.basic_info.n_residues = n_residues
    if (self.basic_info.n_bases == None):
      self.basic_info.n_bases = n_bases


    ## apply resolution cut
    print >> out
    print >> out, "Applying resolution cut"
    print >> out, "-----------------------"

    if self.options.low_resolution is None:
      if self.options.high_resolution is None:
        print >> out, "No resolution cut is made"

    low_cut=float(1e6)
    if self.options.low_resolution is not None:
      low_cut = self.options.low_resolution
      print >> out, "Specified low resolution limit: %3.2f"%(
       float(self.options.low_resolution) )

    high_cut = 0
    if self.options.high_resolution is not None:
      high_cut = self.options.high_resolution
      print >> out, "Specified high resolution limit: %3.2f"%(
       float(self.options.high_resolution) )

    ## perform outlier analyses
    ##
    ## Do a simple outlier analyses please
    print >> out
    print >> out, "Wilson statistics based outlier analyses"
    print >> out, "----------------------------------------"
    print >> out
    native_outlier = data_statistics.possible_outliers(
      miller_array = self.x1,
      prob_cut_ex = self.options.outlier_level_extreme,
      prob_cut_wil = self.options.outlier_level_wilson )
    native_outlier.show(out=out)

    self.x1 = native_outlier.remove_outliers(
      self.x1 )

    ## apply anisotropic scaling  (final B-value will be set to b_add)!
    if self.options.aniso_correction:

      b_final = self.options.b_add
      if b_final is None:
        b_final = 0.0

      print >> out
      print >> out, "Anisotropic absolute scaling of data"
      print >> out, "--------------------------------------"
      print >> out

      aniso_correct = absolute_scaling.ml_aniso_absolute_scaling(
        miller_array = self.x1,
        n_residues = n_residues*\
        self.x1.space_group().order_z()*n_copies_solc,
        n_bases = n_bases*\
        self.x1.space_group().order_z()*n_copies_solc)
      aniso_correct.show(out=out,verbose=1)
      print >> out
      print >> out, "  removing anisotropy for native  "
      print >> out
      u_star_correct_nat = aniso_correct.u_star
      self.x1 = absolute_scaling.anisotropic_correction(
        self.x1,
        aniso_correct.p_scale,
        u_star_correct_nat  )
