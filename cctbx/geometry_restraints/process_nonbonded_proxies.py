from __future__ import division
#from libtbx.utils import Sorry
#from scitbx.array_family import flex
#from libtbx import easy_run
#import iotbx.pdb
#import math
#import sys


class clashes(object):
  """
  Class for clashes
  """
  def __init__(self, clashes_dict):
    """
    clashes_dict  {(iseq, jseq):(distance, sum_vdw_radii)}
    iseq          atom i
    jseq          atom j
    distance      distance between atom i and atom j
    sum_vdw_radii sum of vdW radii
    """
    self.clashes_dict = clashes_dict

  def show(self):
    pass

  def is_clashing(self, iseq):
    pass

  def sort_clashes(self):
    pass



class hbonds(object):
  """
  Class for hbonds
  """
  def __init__(self, hbonds_dict):
    """
    hbonds_dict = {(iseq, jseq, kseq):(H_A_distance, X_A_distance, X_H_A_angle)}
    X-H...A
    iseq          atom X (donor heavy atom)
    jseq          atom H (donor H atom)
    kseq          atom A (acceptor atom)
    H_A_distance
    X_A_distance
    X_H_A_angle
    """
    self.hbonds_dict = hbonds_dict

  def show(self):
    pass

  def forms_hbond(self, iseq):
    pass

  def sort_hbonds(self, sort_distances = True, sort_angles = False):
    pass


class manager():

  def __init__(self,
               model):
    self.model = model
    #

    #
    self._clashes = None
    self._hbonds  = None

    # in manager or do we enfore that input model has H?
    # self._add_H_atoms() ????

  def get_clashes(self):
    """
    Accessor for clashes object
    """
    if not self._clashes:
      self._process_nonbonded_proxies(find_clashes = True)
    else:
      return self._clashes

  def get_hbonds(self):
    """
    Accessor for hbonds object
    """
    if not self._hbonds:
      self._process_nonbonded_proxies(find_hbonds = True)
    else:
      return self._hbonds

  def has_hbonds(self):
    # necessary?
    pass

  def has_clashes(self):
    # necessary?
    pass

  def show(self):
    """
    Print information
    """
    if self.has_clashes():
      self._clashes.show()
    if self.has_hbonds():
      self._hbonds.show()

  def _process_nonbonded_proxies(self,
                                 find_clashes = True,
                                 find_hbonds = False):
    """
    Here is where the calculations are done
    Either all is done at once (clashes, hbonds, other?)
    or it will be modular (use find_clashes and find_hbodns parameters)
    """
    self.geometry_restraints_manager = model.get_restraints_manager().geometry
    self.xrs = model.get_xray_structure()
    self.sites_cart = model.get_sites_cart()
    self.site_labels = xrs.scatterers().extract_labels()
    self.hd_sel = model.get_hd_selection()
    pass
    # instanciate:
    # self._clashes_dict = dict()
    # self._hbonds_dict  = dict()
    # loop, do stuff and fill in the dicts:
    # self._clashes_dict[(iseq, jseq)] = (relevant info)
    # self._hbonds_dict[(iseq, jseq)] = (relevant info)
    # create class:
    # self._clashes = clashes(clashes_dict = clashes_dict)
    # self._hbonds = hbonds(hbonds_dict = hbonds_dict)

