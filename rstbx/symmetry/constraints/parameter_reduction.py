from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints import AGconvert

class symmetrize_reduce_enlarge(object):
                        # symmetrize the metrical matrix &
                        # reduce the number of parameters to reflect symmetry
                        # also, provide back transform to increase
                        # number of parameters back to six

  def __init__(self, space_group):

    self.space_group = space_group
    self.constraints = sgtbx.tensor_rank_2_constraints(
      space_group=self.space_group,reciprocal_space=True)

  def set_orientation(self, orientation, length_unit=1.E-10):
    #provide orientation as either an A matrix (Rossmann) or B matrix (Busing & Levy)
    # data type can be either scitbx.matrix.sqr or scitbx::mat3
    # in either the reciprocal or direct space setting
    # or as a cctbx.crystal_orientation.crystal_orientation
    # if space group is not triclinic the orientation matrix should be close to
    #  symmetrized, but exact symmetrization is done by averaging within the constructor.

    # length unit defaults to 1.E-10 meters = 1 Angstrom

    if "direct_matrix" in dir(orientation):
      self.orientation = orientation # data is already a cctbx orientation
    else:
      from cctbx.crystal_orientation import crystal_orientation
      which_setting = [crystal_orientation(orientation,True),
                       crystal_orientation(orientation,False)]
      #kludgy test for space setting: unit cell volume is never < 40 Angstroms^3
      conversion_to_A3 = (length_unit*length_unit*length_unit)/1.E-30
      select = [a.unit_cell().volume()*conversion_to_A3 > 40.
                for a in which_setting]
      self.orientation = which_setting[select.index(True)]


  def symmetrize(self):
      converter = AGconvert()
      converter.forward(self.orientation) # allows subsequent back-conversion of symmetrized

      # takes member-data orientation; returns symmetrized metrical matrix
      uc = self.orientation.unit_cell()
      avg = self.space_group.average_unit_cell(uc) # assumes direct-space cell
      sym_mm = avg.reciprocal().metrical_matrix()

      converter.validate_and_setG(sym_mm)
      self.orientation = converter.back_as_orientation()

  # it is assumed that metrical_matrix and independent are in reciprocal setting
  def reduce(self,metrical_matrix):
    # takes 6-parameter metrical matrix, returns reduced number of independent parameters
    return self.constraints.independent_params(all_params=metrical_matrix)

  def enlarge(self,independent):
    # takes reduced number independent parameters, returns 6-parameter metrical matrix
    u_star = self.constraints.all_params(independent_params=tuple(independent))
    assert len(u_star) == 6
    return u_star

  def forward_independent_parameters(self):
    # returns the independent parameters given the set_orientation() B matrix
    self.Bconverter=AGconvert()
    self.Bconverter.forward(self.orientation)
    return self.reduce(metrical_matrix = self.Bconverter.G)

  def forward_gradients(self):
    #Specifically for refinement of the B-matrix parameters.
    from rstbx.symmetry.constraints.g_gradients import g_gradients
    gradient_engine = g_gradients(agadaptor = self.Bconverter, symred = self)
    return gradient_engine.dB_dp()

  def backward_orientation(self,independent):
    # given new values of the independent parameters, back-calculate and
    # set the new orientation matrix
    new_mm = self.enlarge(independent)
    self.Bconverter.validate_and_setG(new_mm)
    self.orientation = self.Bconverter.back_as_orientation()
    return self.orientation


if __name__=="__main__":
  from six.moves import cPickle as pickle
  import os
  from labelit.symmetry import metricsym
  from cctbx.sgtbx.bravais_types import bravais_lattice
  from rstbx.symmetry.constraints.g_gradients import finite_difference_test

  this_directory = os.path.dirname(metricsym.__file__)
  F = open(os.path.join(this_directory,"datastore"),"r")
  #datastore represents all subgroup settings in the labelit reference database.

  cases = []
  bravais_types = []
  try:
    while 1:
      cases.append(pickle.load(F))
  except Exception:
    assert len(cases)>900 #should be 980

  for case in cases:
      orient,cs,symmetry = case
      print(orient,cs,symmetry.space_group().type().lookup_symbol())
      S = symmetrize_reduce_enlarge(space_group=symmetry.space_group())
      S.set_orientation(orientation=orient,length_unit=1.E-9)
      S.symmetrize()
      X = S.forward_independent_parameters()
      dB_dp = S.forward_gradients()
      B = S.backward_orientation(independent=X).reciprocal_matrix()
      #continue
      if symmetry.space_group().type().lookup_symbol() not in bravais_types:
        bravais_types.append(symmetry.space_group().type().lookup_symbol())
      finite_difference_test(orient)

  encountered_types = [str(bravais_lattice(i)) for i in bravais_types]
  reference_types = ["aP","mP","mC","oP","oC","oI","oF","tP","tI","hP","hR","cP","cI","cF"]
  encountered_types.sort(); reference_types.sort()
  assert encountered_types == reference_types
  print("OK")
