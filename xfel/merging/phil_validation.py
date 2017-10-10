from __future__ import division
from builtins import object
from libtbx.utils import Sorry

class application(object):
  def __init__(self,param):

    self.param = param
    self.application_level_validation()

  def application_level_validation(self):

    if self.param.merging.reverse_lookup is not None:
      if self.param.data_reindex_op != "h,k,l":
        raise Sorry("The data reindex operator "+self.param.data_reindex_op+
        """ cannot be given in combination with the reverse lookup file
       """+self.param.merging.reverse_lookup+""".  The reverse lookup table itself contains
       a reindexing operator for each image in the data set.""")

class samosa(object):
  def __init__(self,param,allow_model=False):

    self.param = param
    self.allow_model = allow_model
    self.application_level_validation()

  def application_level_validation(self):

    if self.allow_model==False and self.param.model is not None:
      raise Sorry("""For samosa, no PDB structural model is used for frame-to-frame
      scaling, therefore the model phil parameter must be set to None.""")
    if not self.param.scaling.algorithm in ['mark1','levmar']:
      raise Sorry("""Must specify either mark1 or levmar algorithm for scaling.
      (Both algorithms have the same effect within samosa.""")
    #if self.param.significance_filter.apply is True:
    #  raise Sorry("""No significance filter for samosa.  Variance weighting is
    #  used to downweight weak data.""")
    if self.param.raw_data.sdfac_auto or self.param.raw_data.sdfac_refine:
      raise Sorry("""SDFAC manipulation not implemented for samosa.""")
    if self.param.postrefinement.enable:
      raise Sorry("""Postrefinement not currently available in samosa.""")
    if not self.param.include_negatives:
      raise Sorry("""Normally negative values are included in samosa.""")
    if not self.param.scaling.report_ML:
      raise Sorry("""Must have 'report_ML' estimates of mosaic model from integration program""")
    if "PartialityDeff" in self.param.levmar.parameter_flags and \
       "PartialityEtaDeff" in self.param.levmar.parameter_flags:
      raise Sorry("""Only one partiality model""")
