from __future__ import division
from rstbx.array_family import flex
from rstbx.indexing_api import dps_extended
from rstbx.indexing_api.sampling import hemisphere_shortcut
from rstbx.dps_core import Directional_FFT
import math,cmath
from scitbx import matrix
from libtbx.test_utils import approx_equal

class DPS_primitive_lattice:
  def __init__(self, max_cell, recommended_grid_sampling_rad, horizon_phil):
    from libtbx import adopt_init_args
    adopt_init_args(self,locals())

  def set_beam_vector(self,beam):
    self.beam = beam
    self.inv_wave = self.beam.length()

  def set_rotation_axis(self,axis):
    self.axis = axis
    assert self.axis.length() == 1.0

  def set_detector_position(self,origin,d1,d2):
    self.origin = origin
    self.d1 = d1
    self.d2 = d2
    assert self.d1.length() == 1.0
    assert self.d2.length() == 1.0
    assert self.d1.dot(self.d2) == 0.0

  def get_beam_vector_score(self,engine,trial_beam,unique):
    trial_beam = matrix.col(trial_beam)
    nh = min ( engine.getSolutions().size(), 20)
    reciprocal_space_vectors = flex.vec3_double()

    # tile surface to laboratory transformation
    for n in xrange(len(self.raw_spot_input_record)):
      lab_direct = \
      self.origin + self.d1 * self.raw_spot_input_record[n][0] + \
                    self.d2 * self.raw_spot_input_record[n][1]

    # laboratory direct to reciprocal space xyz transformation
      lab_recip = (lab_direct.normalize() * self.inv_wave) - trial_beam

      reciprocal_space_vectors.append ( lab_recip.rotate_around_origin(
        axis=self.axis, angle=self.raw_spot_input_record[n][2], deg=True)
        )

    solutions = engine.getSolutions()
    sum_score = 0.0
    for t in xrange(nh):
      #if t!=unique:continue
      dfft = Directional_FFT(angle = solutions[t], xyzdata = reciprocal_space_vectors,
                            granularity = engine.granularity, amax = engine.amax,
                            F0_cutoff = 11)
      kval = dfft.kval();
      kmax = dfft.kmax();
      kval_cutoff = self.raw_spot_input_record.size()/4.0;
      if ( kval > kval_cutoff ):
        ff=dfft.fft_result;
        projection = trial_beam.dot(matrix.col(solutions[t].dvec))
        kbeam = ((-dfft.pmin)/dfft.delta_p) + 0.5;
        Tkmax = cmath.phase(ff[kmax]);
        backmax = math.cos(Tkmax+(2*math.pi*kmax*kbeam/(2*ff.size()-1)) );
        ### Here it should be possible to calculate a gradient.
        ### Then minimize with respect to two coordinates.  Use lbfgs?  Have second derivatives?
        ### can I do something local to model the cosine wave?
        ### direction of wave travel.  Period. phase.
        sum_score += backmax;
      #if t == unique:
      #  print t, kmax, dfft.pmin, dfft.delta_p, Tkmax,(2*math.pi*kmax*kbeam/(2*ff.size()-1))
    return sum_score

  def index(self,raw_spot_input=None):
    assert raw_spot_input is not None
    self.raw_spot_input_record = raw_spot_input
    # must be x, y, phi in degrees, as a vec3_double

    D = dps_extended()
    D.set_detector(self.d1,self.d2,self.origin)
    D.setMaxcell(self.max_cell)
    reciprocal_space_vectors = flex.vec3_double()

    # tile surface to laboratory transformation
    for n in xrange(len(raw_spot_input)):
      lab_direct = \
      self.origin + self.d1 * raw_spot_input[n][0] + \
                    self.d2 * raw_spot_input[n][1]

    # laboratory direct to reciprocal space xyz transformation
      lab_recip = (lab_direct.normalize() * self.inv_wave) - self.beam

      reciprocal_space_vectors.append ( lab_recip.rotate_around_origin(
        axis=self.axis, angle=raw_spot_input[n][2], deg=True)
        )
      # XXX revisit this.  Why is the angle positive or negative?
    D.setXyzData(reciprocal_space_vectors)
    D.raw_spot_input = raw_spot_input #ersatz until we can refactor this
    D.rotation_vector = self.axis #same comment
    D.wavelength_set = 1./self.inv_wave #same comment
    D.beam_vector = self.beam
    # XXX refactoring strategy:  The DPS_primitive_lattice class renamed to primitive_indexing
    #  then: it inherits directly from dps_extended; is-a rather than has-a

    hemisphere_shortcut(ai = D,
        characteristic_sampling = self.recommended_grid_sampling_rad,
        max_cell = self.max_cell
      )

    scope_check = "Central"
    if scope_check is not None:

      ############  Implement a direct beam check right here #########################
      #score = self.get_beam_vector_score(D,self.beam)
      unique=0
      solutions = D.getSolutions()

      print solutions[unique].dvec
      print self.beam
      # construct two vectors that are perpendicular to the beam.  Gives a basis for refining beam
      beamr0 = self.beam.cross(self.axis).normalize()
      beamr1 = beamr0.cross(self.beam).normalize()
      beamr2 = beamr1.cross(self.beam).normalize()
      print "beamr1",beamr1
      print "beamr2",beamr2
      assert approx_equal(self.beam.dot(beamr1), 0.)
      assert approx_equal(self.beam.dot(beamr2), 0.)
      assert approx_equal(beamr2.dot(beamr1), 0.)
      dv1 = matrix.col(solutions[unique].dvec)
      directionvec = matrix.col((dv1[0],dv1[1],dv1[2]))
      print "directionvec",directionvec
      r1 = beamr1.dot(directionvec)
      r2 = beamr2.dot(directionvec)
      print "r1,r2",r1,r2
      # so the orthonormal vectors are self.beam, beamr1 and beamr2

      grid = 10

      # DO A SIMPLEX MINIMIZATION
      from scitbx.simplex import simplex_opt
      class test_simplex_method(object):
        def __init__(selfOO):
          selfOO.starting_simplex=[]
          selfOO.n = 2
          for ii in range(selfOO.n+1):
            selfOO.starting_simplex.append(flex.random_double(selfOO.n))
          selfOO.optimizer = simplex_opt( dimension=selfOO.n,
                                        matrix  = selfOO.starting_simplex,
                                        evaluator = selfOO,
                                        tolerance=1e-7)
          selfOO.x = selfOO.optimizer.get_solution()


        def target(selfOO, vector):
          newvec = matrix.col(self.beam) + vector[0]*0.0002*beamr1 + vector[1]*0.0002*beamr2
          normal = newvec.normalize() * self.inv_wave
          return -self.get_beam_vector_score(D,normal,unique)

      from scitbx.examples.minimizer_comparisons import derivative
      class test_cma_es(object):
        def __init__(selfOO,l=0):
          selfOO.m = flex.double( [4,4] )
          selfOO.s = flex.double( [2,2])
          selfOO.l = l
          selfOO.fcount = 0
          from cma_es import cma_es_interface
          selfOO.minimizer = cma_es_interface.cma_es_driver( 2, selfOO.m, selfOO.s, selfOO.my_function, selfOO.l )
          print "CMA-ES ITERATIONS", selfOO.minimizer.count, selfOO.fcount,"SOLUTION",  list(selfOO.minimizer.x_final), selfOO.my_function( selfOO.minimizer.x_final )

          selfOO.x = selfOO.minimizer.x_final.deep_copy()


        def compute_functional_and_gradients(selfOO):
          f = selfOO.my_function(selfOO.x)
          g = derivative(selfOO.x, selfOO.my_function,h=1e-4)
          return f, g

        def my_function(selfOO,vector):
          print "CMA-ES",list(vector)
          newvec = matrix.col(self.beam) + vector[0]*0.0002*beamr1 + vector[1]*0.0002*beamr2
          normal = newvec.normalize() * self.inv_wave
          return -self.get_beam_vector_score(D,normal,unique)

      MIN = test_simplex_method()
      #MIN = test_cma_es()
      print "MINIMUM=",list(MIN.x)
      newvec = matrix.col(self.beam) + MIN.x[0]*0.0002*beamr1 + MIN.x[1]*0.0002*beamr2
      normal = newvec.normalize() * self.inv_wave
      self.new_beam = normal
      print "old beam",list(self.beam.elems)
      print "new beam",list(normal.elems)

      scores = flex.double()
      for x in xrange(-grid,grid+1):
       break # skip the grid search
       for y in xrange(-grid,grid+1):
        ref = matrix.col(self.beam)
        newvec = ref + x*0.0002*beamr1 + y*0.0002*beamr2
        normal = newvec.normalize() * self.inv_wave
        scores.append( self.get_beam_vector_score(D,normal,unique) )


      def show_plot(grid,excursi):
        print grid, excursi.size()
        excursi.reshape(flex.grid(grid, grid))

        from matplotlib import pyplot as plt
        plt.figure()
        CS = plt.contour([i*0.2 for i in xrange(grid)],[i*0.2 for i in xrange(grid)], excursi.as_numpy_array())
        plt.clabel(CS, inline=1, fontsize=10, fmt="%6.3f")
        plt.title("Score as to beam likelihood")
        plt.scatter([0.1*(grid-1)],[0.1*(grid-1)],color='g',marker='o')
        #plt.scatter([0.1*(grid-1)+r1],[0.1*(grid-1)+r2],color='g',marker='^')
        plt.scatter([0.1*(grid-1)+0.2*MIN.x[0]] , [0.1*(grid-1)+0.2*MIN.x[1]],color='r',marker='*')
        plt.axes().set_aspect("equal")
        plt.show()

      #show_plot(2 * grid + 1, scores)


    try:
      from rstbx.indexing_api.basis_choice import SelectBasisMetaprocedure as SBM
      pd = {}
      M = SBM(input_index_engine = D,input_dictionary = pd, horizon_phil = self.horizon_phil)

    except Exception,e:
      print e
      import traceback
      traceback.print_exc()
    print "Finished SELECT BASIS with solution M",M
    if 1:#try:

      from rstbx.dps_core.lepage import iotbx_converter
      L = iotbx_converter(D.getOrientation().unit_cell().minimum_cell(),5.0)
      supergroup = L[0]

      triclinic = D.getOrientation().unit_cell()

      cb_op = supergroup['cb_op_inp_best'].c().as_double_array()[0:9]
      orient = D.getOrientation()
      orient_best = orient.change_basis(matrix.sqr(cb_op).transpose())
      constrain_orient = orient_best.constrain(supergroup['system'])
      D.setOrientation(constrain_orient)

      if True:
        for subgroup in L:
          print subgroup.short_digest()
        print "\ntriclinic cell=%s volume(A^3)=%.3f"%(triclinic,triclinic.volume())
        print "\nafter symmetrizing to %s:"%supergroup.reference_lookup_symbol()
        #M.show_rms()
      return D,L
    if 0:#except Exception,e:
      print e
      import traceback
      traceback.print_exc()
