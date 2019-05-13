from __future__ import division
from __future__ import print_function
from six.moves import range
from scitbx.matrix import sqr,col
from math import sin,cos,pi
from scitbx.array_family import flex
import sys
class one_sensor(object):
  def __init__(self,image,sensor,manager):
    self.image = image
    self.sensor = sensor
    self.manager = manager
    self.tiling = self.manager.effective_tiling_as_flex_int(
                    encode_inactive_as_zeroes=True)

    print(list( self.tiling[4*self.sensor[0]:4+4*self.sensor[0]] ))



    grid_radius = 20
    mapp = flex.double(flex.grid(2*grid_radius+1, 2*grid_radius+1))
    print(mapp.focus())

    # not sure if this is the correct beam center; take provisional value
    beam_center = col((float(880.5),float(880.5)))
    gmax = 0.0
    coordmax = (0,0)
    for xi in range(-grid_radius, grid_radius+1):
      for yi in range(-grid_radius, grid_radius+1):
        VV = self.CC(beam_center + col((xi,yi)))
        if VV>gmax:
          gmax = VV
          coordmax = col((xi,yi))
        mapp[(xi+grid_radius,yi+grid_radius)]=VV

    print("max cc %7.4F is at "%gmax, end=' ')
    if False:
      npy = mapp.as_numpy_array()
      from matplotlib import pyplot as plt
      plt.imshow(npy, cmap="hot")
      plt.plot([coordmax[1]+grid_radius],[coordmax[0]+grid_radius],"k.")
      plt.show()

    self.coordmax = coordmax

  def CC(self, beam_center):
    quad = self.tiling[4*self.sensor[0]:4+4*self.sensor[0]]

    asic = self.image.linearintdata.matrix_copy_block(quad[0],quad[1],
                                    quad[2]-quad[0],quad[3]-quad[1])

    #npy = asic.as_numpy_array()
    #from matplotlib import pyplot as plt
    #plt.imshow(npy, cmap="hot")
    #plt.show()

    asci_origin = col((float(quad[0]),float(quad[1])))

    rot45 = sqr((sin(pi/4.),-cos(pi/4.),cos(pi/4.),sin(pi/4.)))

    from xfel.metrology.legacy_scale import quadrant_self_correlation
    min_value = self.image.get_detector()[0].get_trusted_range()[0]
    REF,ROT = quadrant_self_correlation(asic,asci_origin,beam_center,rot45,min_value)
    CCRR = flex.linear_correlation(REF,ROT)

    """initial python implementation
    #rot_asic = flex.double(asic.accessor())
    F0,F1 = asic.focus()

    ref_data = flex.double()
    rot_data = flex.double()
    constant = rot45*(asci_origin - beam_center) +beam_center - asci_origin
    for xcoord in range(quad[2]-quad[0]):
      for ycoord in range(quad[3]-quad[1]):

        acoord = col((float(xcoord),float(ycoord)))
        #prime = rot45*(acoord + asci_origin - beam_center) + beam_center - asci_origin
        prime = rot45*acoord + constant
        prime = (int(round(prime[0],0)),int(round(prime[1],0)))
        if 0<=prime[0]<F0 and 0<=prime[1]<F1:
          #rot_asic[(xcoord,ycoord)]=asic[prime]

          ref_data.append(asic[(xcoord,ycoord)])
          rot_data.append(asic[prime])
    CC = flex.linear_correlation(ref_data,rot_data)

    print "Correlation_coefficient %7.4f %7.4f"%(CC.coefficient(),CCRR.coefficient())
    """

    return CCRR.coefficient()

    #npy = rot_asic.as_numpy_array()
    #from matplotlib import pyplot as plt
    #plt.imshow(npy,cmap="hot")
    #plt.show()

class one_panel(object):
  def __init__(self,image,panel,i_quad,quad,plot=False,multi_angle=True,plot_range=None,show=True):
    self.image = image
    self.panel = panel
    self.i_quad = i_quad
    self.quad = quad
    if multi_angle:
      angles = [20.,22.5, 25.,27.5, 30.,32.5, 35.,37.5,40.,42.5,45.,47.5,50.,52.5,55.,57.5,60.,62.5,65.,67.5,70.]
    else:
      angles = [45.]

    grid_radius = 20
    mapp = flex.double(flex.grid(2*grid_radius+1, 2*grid_radius+1))
    print("Searching a grid with dimensions", mapp.focus())

    beam = image.get_beam()
    beam_center = col(panel.get_beam_centre_lab(beam.get_s0())[0:2])
    gmax = 0.0
    coordmax = (0,0)
    print("Rastering row", end=' ')
    for xi in range(-grid_radius, grid_radius+1):
      print(xi, end=' '); sys.stdout.flush()
      for yi in range(-grid_radius, grid_radius+1):
        delta = col((xi,yi))
        all_VV = flex.double()
        for deg in angles:
          ang_rad = deg * pi/180.
          rotmat = sqr(  (sin(ang_rad), -cos(ang_rad), cos(ang_rad), sin(ang_rad))  )
          all_VV.append(self.CC(beam_center + delta, rotmat))
        #VV = self.CC(beam_center + delta)
        VV = flex.max(all_VV)
        if VV>gmax:
          gmax = VV
          coordmax = delta
        mapp[(xi+grid_radius,yi+grid_radius)]=VV
    print()

    print("max cc %7.4F is at (%d, %d)"%(gmax, coordmax[0], coordmax[1]))
    if plot:
      npy = mapp.as_numpy_array().T # T: transpose
      from matplotlib import pyplot as plt
      fig = plt.figure()
      if plot_range is None:
        vmin = vmax = None
      else:
        vmin, vmax = plot_range
      plt.imshow(npy, cmap="hot", interpolation="nearest", extent=[-grid_radius-1, grid_radius+1, grid_radius+1, -grid_radius-1],
                 vmin = vmin, vmax = vmax)
      ax = plt.colorbar()
      ax.set_label("CC")
      plt.plot([coordmax[0]],[coordmax[1]],"k.")
      plt.title("Rotational autocorrelation of quadrant %d"%i_quad)
      plt.xlabel("X offset (pixels)")
      plt.ylabel("Y offset (pixels)")

      if show:
        plt.show()

    self.coordmax = coordmax
    self.ccmax = gmax

  def CC(self, beam_center, rotmat=None):
    detector = self.image.get_detector()
    angle = [0,3,2,1][self.i_quad] #

    asic = self.image.get_raw_data()[list(detector.get_names()).index(self.panel.get_name())].matrix_rot90(angle)

    p_w, p_h = self.panel.get_image_size()
    b = [self.panel.get_pixel_lab_coord((0    ,0    )),
         self.panel.get_pixel_lab_coord((p_w-1,0    )),
         self.panel.get_pixel_lab_coord((p_w-1,p_h-1)),
         self.panel.get_pixel_lab_coord((0    ,p_h-1))]
    asic_origin = col(self.panel.millimeter_to_pixel((min([p[0] for p in b]),
                                                      min([p[1] for p in b]))))

    if rotmat is None:
      rot45 = sqr((sin(pi/3.),-cos(pi/3.),cos(pi/3.),sin(pi/3.)))
    else: rot45 = rotmat

    from xfel.metrology.legacy_scale import quadrant_self_correlation
    min_value = self.image.get_detector()[0].get_trusted_range()[0]
    REF,ROT = quadrant_self_correlation(asic.as_double(),asic_origin,beam_center,rot45,min_value)
    CCRR = flex.linear_correlation(REF,ROT)

    return CCRR.coefficient()

    #npy = rot_asic.as_numpy_array()
    #from matplotlib import pyplot as plt
    #plt.imshow(npy,cmap="hot")
    #plt.show()
