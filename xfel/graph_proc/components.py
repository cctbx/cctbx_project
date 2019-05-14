""" Module for treating images as the vertices of a graph. Includes both Edge and Vertex (ImageNode) classes. """
from __future__ import absolute_import, division, print_function
from scitbx.matrix import sqr
from cctbx.uctbx import unit_cell
import numpy as np
import math
import logging
from cctbx.array_family import flex
from xfel.clustering.singleframe import SingleFrame

class ImageNode(SingleFrame):
  """ Extends a SingleFrame object to make it a vertex in a graph. Adds a list
   of edges.
  """

  def __init__(self, *args, **kwargs):
    """ Constructor is same as for SingleFrame object, but has additional
    kwargs:
    :param kwargs['scale']: default True. Specifies if the images should be scaled upon creation. Mainly switched off for testing.
    :param kwargs['use_b']: default True. If false, only initialise the scale factor, not the B factor.
    """
    SingleFrame.__init__(self, *args, **kwargs)
    if hasattr(self, 'miller_array'):  # i.e. if the above worked.
      self.use_scales = kwargs.get('scale', True)
      self.use_b = kwargs.get('use_b', True)
      self.edges = []  # this is populated when the Graph is made.
      self.partialities = self.calc_partiality(self.get_x0(),
                                               update_wilson=self.use_scales)
      self.scales = self.calc_scales(self.get_x0())
      self.label = None  # to be used for classification after instantiation
      self.source = None  # for testing, when we know the 'source' of the image
      self.params = self.get_x0()


  def calc_partiality(self, params_in, update_wilson=False):
    """ Reuturn a partiality vector for all the reflections in self, using the
    parameters defined in the array params.

    :params: a tuple of the form appropriate for the crystal symetry, such as
    one produced by get_x0()
    :return: a list of partialities for the miller indicies in self.
    """
    from xfel.cxi.postrefine.mod_leastsqr import prep_output, \
      get_crystal_orientation
    from xfel.cxi.postrefine.mod_partiality \
      import partiality_handler as PartialityHandler

    # Expand to triclinic params, and unpack
    params_all = prep_output(params_in, self.crystal_system)
    G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = params_all

    # Create a unit cell based on params, not the initialised values.
    try:
      uc = unit_cell((a, b, c, alpha, beta, gamma))
    except RuntimeError:
      return  None
      logging.warning("Could not create unit cell with ({}, {}, {}, {}, {}, {})"
                      ". Resetting unit cell to original value."
                      .format(a, b, c, alpha, beta, gamma))
      uc = self.orientation.unit_cell()

    ortho_matrix = uc.orthogonalization_matrix()
    orientation_matrix = self.orientation.crystal_rotation_matrix()

    crystal_init_orientation = get_crystal_orientation(ortho_matrix,
                                                       orientation_matrix)
    crystal_orientation_model = crystal_init_orientation \
      .rotate_thru((1, 0, 0), rotx).rotate_thru((0, 1, 0), roty)
    a_star = sqr(crystal_orientation_model.reciprocal_matrix())

    # set up variables
    miller_indices = self.miller_array.indices()
    ph = PartialityHandler(self.wavelength, None)
    bragg_angle_set = self.miller_array.two_theta(wavelength=self.wavelength) \
      .data()
    mm_predictions = self.pixel_size * self.mapped_predictions
    alpha_angle_obs = flex.double([math.atan(abs(pred[0] - self.xbeam) /
                                             abs(pred[1] - self.ybeam))
                                   for pred in mm_predictions])
    assert len(alpha_angle_obs) == len(self.miller_array.indices()), \
      'Size of alpha angles and observations are not equal %6.0f, %6.0f' \
      % (len(alpha_angle_obs), len(self.miller_array.indices()))

    partiality = ph.calc_partiality_anisotropy_set(a_star, miller_indices, ry,
                                                   rz, re, bragg_angle_set,
                                                   alpha_angle_obs)
    if update_wilson:
      self.minus_2B, self.G, self.log_i, self.sinsqtheta_over_lambda_sq, \
      self.wilson_err = self.init_calc_wilson(use_b_factor=self.use_b,
                                              i_corrections=(1 / partiality))
    return partiality

  def calc_scales(self, params_in):
    """ Calculate an array of scales based on scale, B and wavelength using the
    equation $scale * exp(-2*B*(sin(theta)/wavelength)^2)$

    Reuturn a scale vector for all the reflections in self, using the
    parameters defined in the array params.

    :params: a tuple of the form appropriate for the crystal symetry, such as
    one produced by get_x0(). This method only uses params[0] (scale) and
    params[1] (B)

    :return: a list of scales for all the miller indicies in self
    """
    if self.use_scales:
      scale = params_in[0]
      B = params_in[1]
      sin_sq_theta = self.miller_array.two_theta(wavelength=self.wavelength) \
        .sin_theta_over_lambda_sq().data()

      scales = scale * self.miller_array.data()
      exp_arg = flex.double(-2 * B * sin_sq_theta)
      return flex.double(flex.double(scales) * flex.exp(exp_arg))
    else:
      # Horrible way to get vector of ones...
      return flex.double(self.miller_array.data()/self.miller_array.data())


  def get_x0(self):
    """
    Return the initial estimages of parameters to be refined for this frame
    :return: inital (G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma)
    for lowest symetry, up to (G, B, rotx, roty, ry, rz, re, a) for cubic.
    values, using:

      - G: initial guess from linear regression
      - B: initial guess from linear regression
      - rotx: 0
      - roty: 0
      - ry: spot radius initial guess
      - rz: spot radius initial guess
      - re: 0.0026 (From Winkler et al. paper)
      - [a, b, c, alpha, beta, gamma]: from initial indexing solution, dependent on crystal system:

            Triclinic: (G, B, rotx, roty, ry, rz, re, a)
            Monoclinic: (G, B, rotx, roty, ry, rz, re,a,b,c,beta])
            Orthorhombic: (G, B, rotx, roty, ry, rz, re,a,b,c)
            Tetragonal: (G, B, rotx, roty, ry, rz, re,a,c)
            Trigonal or Hexagonal: (G, B, rotx, roty, ry, rz, re,a,c)
            Cubic: (G, B, rotx, roty, ry, rz, re,a)
    """
    from xfel.cxi.postrefine.mod_leastsqr import prep_input
    from xfel.cxi.postrefine.test_rs import calc_spot_radius

    a_star = sqr(self.orientation.reciprocal_matrix())
    miller_indices = self.miller_array.indices()
    spot_radius = calc_spot_radius(a_star, miller_indices, self.wavelength)
    x_init = [self.G, - 1 * self.minus_2B / 2, 0, 0,
              spot_radius, spot_radius, 0.0026]
    x_init.extend(self.uc)
    x0_all = np.array(x_init)
    x0 = prep_input(x0_all, self.crystal_system)
    return x0


class Edge:
  """
  .. note::
    Developmental code. Do not use without contacting zeldin@stanford.edu

  Defines an undirected edge in a graph. Contains the connecting vertices, and a weight.
  """

  def __init__(self, vertex_a, vertex_b, weight):
    self.vertex_a = vertex_a
    self.vertex_b = vertex_b
    self.weight = weight
    self.intra = None  # True means we assume it connects two edges from the
                       # same class, False means it is an iter-class edge.

  def mean_residual(self):
    """
    :return: the mean residual of the absolute value of the all the weigts on this edge.
    """
    return np.mean(np.abs(self.residuals()))

  def other_vertex(self, vertex):
    """
    Simple method to get the other vertex along and edge.

    :param vertex: a vertex that is on one end of this edge
    :return: the vertex at the other end of the edge
    """
    assert vertex == self.vertex_a or vertex == self.vertex_b
    if vertex is self.vertex_a:
      return self.vertex_b
    elif vertex is self.vertex_b:
      return self.vertex_a

  @staticmethod
  def _calc_residuals(va, vb, pa, pb, sa, sb):

    mtch_indcs = va.miller_array.match_indices(vb.miller_array,
                                               assert_is_similar_symmetry=False)

    va_selection = mtch_indcs.pair_selection(0)
    vb_selection = mtch_indcs.pair_selection(1)

    sp_a = pa.select(va_selection) * sa.select(va_selection)
    sp_b = pb.select(vb_selection) * sb.select(vb_selection)

    ia_over_ib = va.miller_array.data().select(va_selection) / \
                 vb.miller_array.data().select(vb_selection)

    residuals = (flex.log(sp_a) - flex.log(sp_b) - flex.log(ia_over_ib))
    residuals = residuals.as_numpy_array()
    #logging.debug("Mean Residual: {}".format(np.mean(residuals)))
    return residuals[~np.isnan(residuals)]

  def residuals(self):
    """
    Calculates the edge residual, as defined as the sum over all common miller indices of:
      log(scale * partiality of a) - log(scale * partiality of b) - log(I_a/I_b)

    :return: the residual score for this edge
    """
    # 1. Create flex selection array
    # 2. Trim these so that they only contain common reflections
    # 3. Calculate residual

    partialities_a = self.vertex_a.partialities
    partialities_b = self.vertex_b.partialities
    scales_a = self.vertex_a.scales
    scales_b = self.vertex_b.scales

    return Edge._calc_residuals(self.vertex_a, self.vertex_b,
                           partialities_a, partialities_b,
                           scales_a, scales_b)
