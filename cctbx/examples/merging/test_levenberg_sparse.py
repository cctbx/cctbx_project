from __future__ import division
import os
import math
from libtbx.development.timers import Timer
from scitbx.array_family import flex
from cctbx.examples.merging.task4 import prepare_simulation_with_noise
from cctbx.examples.merging.data_utilities import I_and_G_base_estimate, plot_it, show_correlation
from cctbx.examples.merging.data_subset import mapper_factory
from scitbx.lstbx import normal_eqns
from scitbx.lstbx import normal_eqns_solving
from scitbx.matrix import sqr,col

from cctbx.uctbx import unit_cell
uc_params = (281,281,165,90,90,120)
uc = unit_cell(uc_params)

from cctbx.examples.merging import xscale6e as base_class
class levenberg_helper(base_class, normal_eqns.non_linear_ls_mixin):
  def __init__(pfh, initial_estimates):
    super(levenberg_helper, pfh).__init__(n_parameters=len(initial_estimates))
    pfh.x_0 = flex.double(initial_estimates)
    pfh.restart()
    pfh.counter = 0
    # do this outside the factory function pfh.set_cpp_data(fsim, NI, NG)

  def restart(pfh):
    pfh.x = pfh.x_0.deep_copy()
    pfh.old_x = None

  def step_forward(pfh):
    pfh.old_x = pfh.x.deep_copy()
    pfh.x += pfh.step()

  def step_backward(pfh):
    assert pfh.old_x is not None
    pfh.x, pfh.old_x = pfh.old_x, None

  def parameter_vector_norm(pfh):
    return pfh.x.norm()

  def build_up(pfh, objective_only=False):
    if not objective_only: pfh.counter+=1
    pfh.reset()
    pfh.access_cpp_build_up_directly_eigen_eqn(objective_only, current_values = pfh.x)

class xscale6e(object):

  def __init__(self,Ibase,Gbase,FSIM,curvatures=False):
    self.counter = 0
    self.x = flex.double(list(Ibase) + list(Gbase) + list(flex.double(len(Gbase), 0.0)))
    self.N_I = len(Ibase)
    self.N_G = len(Gbase)
    self.N_raw_obs = FSIM.raw_obs.size()
    print "# structure factors:",self.N_I, "# frames:",self.N_G

    self.helper = levenberg_helper(initial_estimates = self.x)
    self.helper.set_cpp_data(FSIM, self.N_I, self.N_G)
    self.helper.restart()

    iterations = normal_eqns_solving.levenberg_marquardt_iterations_encapsulated_eqns(
               non_linear_ls = self.helper,
               n_max_iterations = 5000,
               track_all=True,
               step_threshold = 0.0001
    )
    print "End of minimisation: Converged", self.helper.counter,"cycles"
    chi_squared = self.helper.objective() * 2.
    print "obj",chi_squared
    print "# of obs:",FSIM.raw_obs.size()
    dof = FSIM.raw_obs.size() - ( self.N_I + 2 * self.N_G )
    print "degrees of freedom =",dof
    print "chisq/dof: %7.3f"%(chi_squared / dof)
    print

  def packed_to_all(self,packed):
    Nx = len(self.helper.x)
    all_elems = flex.double(Nx*Nx)
    ctr = 0
    for x in xrange(Nx):
      x_0 = ctr
      for y in xrange(Nx-1,x-1,-1):
        all_elems[ Nx*x+y ] = packed[x_0+(y-x)]
        ctr += 1
        if x!= y:
          all_elems[ Nx*y+x ] = packed[x_0+(y-x)]
    return all_elems
  def pretty(self, matrix, max_col=67,format="%6.0g",zformat="%6.0f"):
    Nx = len(self.helper.x)
    for islow in xrange(Nx):
      for ifast in xrange(min(max_col,Nx)):
        if matrix(islow,ifast)==0.:
          print zformat%0,
        else:
          print format%(matrix(islow,ifast)),
      print
    print
  def prettynz(self, matrix, max_col=67,format="%6.0g",zformat="%6.0f"):
    Nx = len(self.helper.x)
    for x in xrange(Nx):
      for y in xrange(min(max_col,Nx)):
          if abs(matrix(x,y))< 1E-10:
            print zformat%0,
          else:
            print format%(matrix(x,y)),
      print
  def permutation_ordering_to_matrix(self,ordering):
    matcode = flex.int(ordering.size()**2)
    Nx = ordering.size()
    for islow in xrange(Nx):
      matcode[Nx*islow + ordering[islow]] = 1
    return sqr(matcode)
  def lower_triangular_packed_to_matrix(self,lower):
    Nx = len(self.helper.x)
    lower_elem = flex.double(Nx*Nx)
    ctr = 0
    for irow in xrange(Nx): # loop over rows
      for icol in xrange(irow + 1): # loop over columns
        lower_elem[Nx*irow + icol] = lower[ctr]
        ctr+=1
    return sqr(lower_elem)
  def diagonal_vector_to_matrix(self,diag):
    diag_elem = flex.double(diag.size()**2)
    for irow in xrange(diag.size()): # loop over rows
      diag_elem[diag.size()*irow + irow] = diag[irow]
    return sqr(diag_elem)
  def unstable_matrix_inversion_diagonal(self,Lower,Diag,Transpose):
    ### Can't use the Cholesky factorization to derive the Variance-Covariance matrix
    ### Demonstrate that inverting the matrix is numerically unstable
    Nx = len(self.helper.x)
    error_diagonal_elems = flex.double(Nx)
    for j_element in xrange(Nx):

      # now solve for the vector p = D * LT * x by substitution in eqn L * p = b
      # b is the column vector of zeroes except jth_element is 1.
      p = flex.double(Nx)
      p[j_element] = 1.
      for p_idx in xrange(j_element+1,Nx):
        for subs_idx in xrange(j_element,p_idx):
          p[p_idx] -= Lower(p_idx,subs_idx) * p[subs_idx]

      Pvec = col(p)

      # now solve for the vector q = LT * x by division in eqn D * q = b
      q = flex.double([ p[i] / Diag(i,i) for i in xrange(Nx)] )
      #  this is the unstable step.  We can't divide by tiny denominators
      Qvec = col(q)

      # now solve for the jth element of x in the eqn LT * x = q
      xelem = flex.double(Qvec.elems)

      for x_idx in xrange(Nx-1,j_element-1,-1):  #comment this in for production
      #for x_idx in xrange(Nx-1,-1,-1):
        for subs_idx in xrange(x_idx+1, Nx):
          xelem[x_idx] -= Transpose(x_idx,subs_idx) * xelem[subs_idx]
      Xvec = col(xelem) # got the whole vector; only need j_element for the error matrix diagonal
      error_diagonal_elems[j_element] = xelem[j_element]
    return col(error_diagonal_elems)

  def esd_plot(self):
    print "OK esd"

    ### working on the esd problem:
    self.helper.build_up()
    norm_mat_packed_upper = self.helper.get_normal_matrix()
    norm_mat_all_elems = self.packed_to_all(norm_mat_packed_upper)
    diagonal_curvatures = self.helper.get_normal_matrix_diagonal()
    NM = sqr(norm_mat_all_elems)

    print "The normal matrix is:"
    self.pretty(NM)

    from scitbx.linalg.svd import inverse_via_svd
    svd_inverse,sigma = inverse_via_svd(NM.as_flex_double_matrix())

    print "ia",len(svd_inverse),len(sigma)
    IA = sqr(svd_inverse)
    for i in xrange(self.helper.x.size()):
      if i == self.N_I or i == self.N_I + self.N_G:  print
      print "%2d %10.4f %10.4f %10.4f"%(
        i, self.helper.x[i], math.sqrt(1./diagonal_curvatures[i]), math.sqrt(IA(i,i)))

    from matplotlib import pyplot as plt
    plt.plot(flex.sqrt(flex.double([IA(i,i) for i in xrange(self.N_I)])),
             flex.sqrt(1./diagonal_curvatures[:self.N_I]), "r.")
    plt.title("Structure factor e.s.d's from normal matrix curvatures vs. SVD variance diagonal")
    plt.axes().set_aspect("equal")
    plt.show()
    return

    # additional work to validate the Cholesky factorization and investigate stability:
    identity = IA * NM

    print "verify identity:"
    self.pretty(identity,max_col=58,format="%7.1g")
    # we can fool ourselves that the SVD gave us a perfect inverse:
    self.pretty(identity,max_col=72,format="%4.0f")

    ### figure out stuff about permutation matrices
    self.helper.build_up()

    ordering = self.helper.get_eigen_permutation_ordering()
    print "ordering:",list(ordering)
    matcode = self.permutation_ordering_to_matrix(ordering)

    print "matcode:"
    self.pretty(matcode,max_col=72,format="%1d",zformat="%1d")

    permuted_normal_matrix = (matcode.inverse())* NM *matcode
    print "product"
    self.pretty(permuted_normal_matrix)

    ### Now work with the Cholesky factorization
    cholesky_fac_packed_lower = self.helper.get_cholesky_lower()
    Lower = self.lower_triangular_packed_to_matrix(cholesky_fac_packed_lower)
    print "lower:"
    self.pretty(Lower,max_col=59,format="%7.0g",zformat="%7.0g")

    Transpose = Lower.transpose()
    print "transpose"
    self.pretty(Transpose,max_col=59,format="%7.0g",zformat="%7.0g")

    diagonal_factor = self.helper.get_cholesky_diagonal()
    Diag = self.diagonal_vector_to_matrix(diagonal_factor)
    print "diagonal:"
    self.pretty(Diag,max_col=59,format="%7.0g",zformat="%7.0g")

    Composed = Lower * Diag * Transpose
    print "composed"
    self.pretty(Composed,max_col=67,format="%6.0g",zformat="%6.0g")

    Diff = Composed - permuted_normal_matrix
    print "diff"
    self.prettynz(Diff,max_col=67,format="%6.0g",zformat="%6.0g")
    #  OK, this proves that L * D * LT = P * A * P-1
    #  in other words, Eigen has correctly factored the permuted normal matrix
    ############

    Variance_diagonal = self.unstable_matrix_inversion_diagonal(Lower,Diag,Transpose)

    for i in xrange(self.helper.x.size()):
      if i == self.N_I or i == self.N_I + self.N_G:  print
      print "%2d %10.4f %10.4f %10.4f"%(
        i, self.helper.x[i], math.sqrt(1./diagonal_curvatures[i]), math.sqrt(IA(i,i))),
      print "svd err diag: %10.4f"%(IA(i,i)),"eigen: %15.4f"%(Variance_diagonal[ordering[i]])


  def unpack(self):
    return (
      self.helper.x[ 0 : self.N_I ],                  # fitted intensity
      self.helper.x[ self.N_I : self.N_I + self.N_G ],# fitted G
      self.helper.x[ self.N_I + self.N_G : ]          # fitted B
    )

  def unpack_stddev(self):
    # the data-to_parameter ratio will control which method for returning e.s.d's
    data_to_parameter = float(self.N_raw_obs) / self.helper.x.size()
    self.helper.build_up()
    if data_to_parameter <= 4. and self.helper.x.size() < 500:
      # estimate standard deviations by singular value decomposition
      norm_mat_packed_upper = self.helper.get_normal_matrix()
      norm_mat_all_elems = self.packed_to_all(norm_mat_packed_upper)
      NM = sqr(norm_mat_all_elems)
      from scitbx.linalg.svd import inverse_via_svd
      svd_inverse,sigma = inverse_via_svd(NM.as_flex_double_matrix())
      IA = sqr(svd_inverse)
      estimated_stddev = flex.double([math.sqrt(IA(i,i)) for i in xrange(self.helper.x.size())])
    else:
      # estimate standard deviations by normal matrix curvatures
      diagonal_curvatures = self.helper.get_normal_matrix_diagonal()
      estimated_stddev = flex.sqrt(1./diagonal_curvatures)
    return (
      estimated_stddev[ 0 : self.N_I ],                  # fitted intensity
      estimated_stddev[ self.N_I : self.N_I + self.N_G ],# fitted G
      estimated_stddev[ self.N_I + self.N_G : ]          # fitted B
    )

  def show_summary(self):
    print "%d cycles"%self.counter
    self.helper.show_eigen_summary()

class execute_case(object):
 def __init__(self,datadir,n_frame,transmittance,apply_noise,plot=False,esd_plot=False,half_data_flag=0):
  # read the ground truth values back in
  import cPickle as pickle
  ordered_intensities = pickle.load(open(os.path.join(datadir,"intensities.pickle"),"rb"))
  frames = pickle.load(open(os.path.join(datadir,"frames.pickle"),"rb"))

  sim = pickle.load(open(os.path.join(datadir,"simulated%05d_0.pickle"%n_frame),"rb"))
  print "accepted obs %d"%(len(sim["observed_intensity"]))

  FSIM = prepare_simulation_with_noise(sim, transmittance=transmittance,
                                       apply_noise=apply_noise,
                                       ordered_intensities=ordered_intensities,
                                       half_data_flag=half_data_flag)

  I,I_visited,G,G_visited = I_and_G_base_estimate(FSIM)
  model_I = ordered_intensities.data()[0:len(I)]
  model_G = frames["scale_factors"][0:len(G)]
  model_B = frames["B_factors"][0:len(G)]

  T = Timer("%d frames, %f transmittance, %s noise"%(
             n_frame, transmittance, {False:"NO", True:"YES"}[apply_noise]))

  mapper = mapper_factory(xscale6e)
  minimizer = mapper(I,G,I_visited,G_visited,FSIM)

  del T
  minimizer.show_summary()

  Fit_I, Fit_G, Fit_B = minimizer.e_unpack()
  show_correlation(Fit_G,model_G,G_visited,"Correlation of G:")
  show_correlation(Fit_B,model_B,G_visited,"Correlation of B:")
  show_correlation(Fit_I,model_I,I_visited,"Correlation of I:")
  Fit_I_stddev, Fit_G_stddev, Fit_B_stddev = minimizer.e_unpack_stddev()

  if plot:
    plot_it(Fit_G, model_G)
    plot_it(Fit_B, model_B,mode="B")
    plot_it(Fit_I, model_I, mode="I")
  print

  if esd_plot:
    minimizer.esd_plot()

  from cctbx.examples.merging.show_results import show_overall_observations
  table1,self.n_bins,self.d_min = show_overall_observations(
           Fit_I,Fit_I_stddev,model_I,I_visited,
           ordered_intensities,frames,FSIM,title="Statistics for all reflections")

  self.FSIM=FSIM
  self.ordered_intensities=ordered_intensities
  self.Fit_I=Fit_I
  self.Fit_I_stddev=Fit_I_stddev
  self.I_visited=I_visited

if __name__=="__main__":

  datadir = os.path.join(os.environ["HOME"],"rosie_xds","xscale_reserve") # Get files directly from author, NKS
  plot_flag=False
  execute_case(datadir, n_frame=400, transmittance=0.00001, apply_noise=True, plot=plot_flag)
  print "OK"
