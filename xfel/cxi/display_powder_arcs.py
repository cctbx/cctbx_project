from __future__ import absolute_import, division, print_function
from six.moves import range
import math
import mmtbx.programs.fmodel
import mmtbx.utils
import iotbx.pdb.fetch
from iotbx import pdb
from cctbx.array_family import flex


def is_pinwheel_region(dx,dy):
  angledeg = math.atan2(dy,dx)*180./math.pi
  return (30.<angledeg<60.) or (120.<angledeg<150.) or \
         (-30.>angledeg>-60.) or (-120.>angledeg>-150.)

def get_mmtbx_icalc(code,d_min, anomalous_flag = False):

    #pdbtools.run([file, "high_resolution=%f"%self.d_min, "label=f-obs-nks"
    #              "k_sol=0.35", "b_sol=60", "b_cart=1 2 -3 0 0 0", "--f_model", "r_free=0.1"])
    # get xray_structure

    pdb_url = iotbx.pdb.fetch.fetch(id=code)
    pdb_input = pdb.input(source_info=None,lines=pdb_url.readlines())

    xray_structure = pdb_input.xray_structure_simple()

    phil2 = mmtbx.programs.fmodel.master_phil
    params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    ISO_ALLOWANCE = 0.1 # isomorphous recip cell volume changes no more than 10%
    params2.high_resolution = d_min / math.pow( (1.+ISO_ALLOWANCE),(1./3.) )
    params2.output.type = "real"
    if True :
      params2.fmodel.k_sol = 0.35
      params2.fmodel.b_sol = 46
    f_model = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = xray_structure,
      f_obs          = None,
      add_sigmas     = True,
      params         = params2).f_model
    i_model = f_model.as_intensity_array()


    return  i_model.map_to_asu()

def apply_gaussian_noise(image,params):
    from scitbx.random import variate,normal_distribution
    import numpy
    G = variate(normal_distribution(mean=2.0,sigma=0.5))
    gaussian_noise = flex.double(G(image.linearintdata.size()))
    #image.linearintdata += flex.int(gaussian_noise.as_numpy_array().astype(numpy.int32))
    image.linearintdata += gaussian_noise

def superimpose_powder_arcs(image,params):

    #put in some random noise to improve display
    image.read()
    apply_gaussian_noise(image,params)

    pxlsz = image.pixel_size
    distance = image.distance
    wavelength = image.wavelength
    eps = 0.01
    beamx = eps+image.beamx/pxlsz
    beamy = eps+image.beamy/pxlsz
    print(beamx,beamy)

    image.read()
    data = image.linearintdata
    detector_d = flex.double()
    detector_radius = flex.double()

    for x in range(data.focus()[0]):
      dx = x-beamx
      dx_sq = dx*dx;

      for y in range(data.focus()[1]):
        dy = y-beamy
        radius = math.sqrt(dx_sq+dy*dy)*pxlsz
        detector_radius.append(radius)
        theta = 0.5 * math.atan(radius/distance)
        resol_d = wavelength/(2.*math.sin(theta))

        if is_pinwheel_region(dx,dy):
          detector_d.append(resol_d)
        else:
          detector_d.append(0.0)
    print(detector_d.size())
    d_order_image = flex.sort_permutation(detector_d,reverse=True)

    #code = "2oh5"
    code = params.viewer.powder_arcs.code
    d_min = 2.0
    sf = get_mmtbx_icalc(code,d_min)
    sf.show_summary()

    cell = sf.unit_cell()
    millers = sf.indices()
    FUDGE_FACTOR = 0.02
    intensities = FUDGE_FACTOR*sf.data()
    spot_d = cell.d(millers)
    spot_two_theta = cell.two_theta(millers,wavelength=wavelength)
    spot_radius = distance * flex.tan(spot_two_theta)

    d_order_spots = flex.sort_permutation(spot_d,reverse=True)
    print(list(d_order_spots))

    spot_ptr_min = 0
    spot_ptr_max = 1

    THREESIGMA=1.5
    SIGMA = 0.20
    SCALE = 0.003
    for x in range(50000):#len(detector_d)):
      if x%10000==0: print(x)
      if detector_d[d_order_image[x]]>2.1:
        pxl_radius = detector_radius[d_order_image[x]]
        #print "pixel radius",pxl_radius,"resolution",detector_d[d_order_image[x]]
        if True:# data[d_order_image[x]]>0:
          data[d_order_image[x]]=0
          # increase spot_ptr_max???
          while spot_radius[d_order_spots[spot_ptr_max]]-THREESIGMA < pxl_radius:
            spot_ptr_max += 1
            print(x, pxl_radius,"increase spot_ptr_max" , spot_ptr_max, spot_radius[d_order_spots[spot_ptr_max]]-THREESIGMA)
          while spot_radius[d_order_spots[spot_ptr_min]]+THREESIGMA < pxl_radius:
            spot_ptr_min += 1
            print(x, pxl_radius,"increase spot_ptr_min" , spot_ptr_min)
          for sidx in range(spot_ptr_min,spot_ptr_max):
            spot_rad = spot_radius[d_order_spots[sidx]]
            delta_rad = spot_rad-pxl_radius
            #gaussian distribution with SIGMA=1 (1 pixel sigma)
            #divide through by spot_radius since energy has to be spread around increasing ring
            #  circumference
            #print "changing data ",x,d_order_image[x]
            data[d_order_image[x]] += int(SCALE*(1./(math.sqrt(2.*math.pi*SIGMA*SIGMA)))*math.exp(
              -0.5*delta_rad*delta_rad/(SIGMA*SIGMA))*intensities[d_order_spots[sidx]]/spot_rad)
