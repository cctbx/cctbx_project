from __future__ import division
from libtbx.utils import Sorry

class application:
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
