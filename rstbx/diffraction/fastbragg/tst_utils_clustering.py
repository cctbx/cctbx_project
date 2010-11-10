from math import sqrt,pi,acos,sin,atan
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx.matrix import col,sqr

def specific_libann_cluster(data,intensity_cutoff = 25,distance_cutoff=17):
  from annlib_ext import AnnAdaptor
  #construct a new data structure containing only pixels > intensity_cutoff
  #input data structure is a flex array
  info = {}
  shape = data.accessor().focus()
  for slow in xrange(shape[0]):
    for fast in xrange(shape[1]):
      if data[(slow,fast)] > intensity_cutoff:
        info[(slow,fast)] = data[(slow,fast)]

  Ktree = int(distance_cutoff*distance_cutoff*pi)
  Sq_cut = distance_cutoff*distance_cutoff # distance < distance_cutoff => the two points are clustered

  all_pixels = flex.double()
  all_keys = info.keys()
  for key in info.keys():
    all_pixels.append(key[0]); all_pixels.append(key[1])

  distance_tree = AnnAdaptor(data=all_pixels,dim=2,k=Ktree)
  distance_tree.query(all_pixels)
  clusters = []
  membership_lookup = {}

  for i_query_pt in xrange(len(all_keys)):
      query_coords = all_keys[i_query_pt]
      query_ok = True
      for i_target_pt in xrange(Ktree):
        target_coords = all_keys[distance_tree.nn[Ktree*i_query_pt+i_target_pt]]
        if distance_tree.distances[Ktree*i_query_pt+i_target_pt] < Sq_cut:
          if info[query_coords] < info[target_coords]:
            query_ok = False
            break
      if query_ok:
        membership_lookup[query_coords]=info[query_coords]

  return membership_lookup

def index_wrapper(positions,info,pdb_object):
  from rstbx.dps_core import dps_core
  from rstbx.dps_core.sampling import HemisphereSamplerBase as HemisphereSampler
  from cctbx.uctbx import unit_cell

  uc = pdb_object.unit_cell()

  sampling = get_recommended_sampling(pdb_object, info, positions)

  raw_spot_input = flex.vec3_double()
  pixel_sz = info.D.pixel_sz
  for pos in positions:
    raw_spot_input.append((pos[0]*pixel_sz, pos[1]*pixel_sz, 0.0))

  #convert raw film to camera, using labelit coordinate convention
  camdata = flex.vec3_double()
  auxbeam = col((info.C.Ybeam,info.C.Zbeam,0.0));
  film_2_camera = sqr((-1,0,0,0,-1,0,0,0,1)).inverse();
  for x in xrange(len(raw_spot_input)):
    camdata.append( auxbeam + film_2_camera * col(raw_spot_input[x]) )

  #convert camera to reciprocal space xyz coordinates
  xyzdata = flex.vec3_double()

  for x in xrange(len(camdata)):
    cam = col(camdata[x])
    auxpoint = col((cam[0],cam[1],info.C.distance));
    xyz = ( auxpoint / (info.C.lambda0*1E10 * auxpoint.length()) );
    xyz = xyz - col((0.0, 0.0, 1.0/(info.C.lambda0*1E10))); #translate recip. origin

    xyzdata.append( xyz );

  core_ai = dps_core()
  core_ai.setXyzData(xyzdata)
  core_ai.setMaxcell(1.25*max(uc.parameters()[0:3]))
  """
add select files to svn status
check this in to rstbx!!
email Holton
email AZET"""

  H = HemisphereSampler(
      characteristic_grid = sampling,
      max_cell=1.25*max(uc.parameters()[0:3]))
  H.hemisphere(core_ai,size=30,cutoff_divisor=4.) # never change these parameters

  from rstbx.dps_core.basis_choice import SelectBasisMetaprocedure as SBM
  M = SBM(core_ai)

  return core_ai,uc

def get_recommended_sampling(pdb_obj,info,spots):
    p1_uc = pdb_obj.unit_cell()
    largest_p1_cell_dimension = max(p1_uc.parameters()[0:3])

    relevant_resolution =  info.C.lambda0*1E10

    pixel_sz = info.D.pixel_sz
    xtd=info.C.distance
    xbeam=info.C.Ybeam
    ybeam=info.C.Zbeam
    resolutions=[]
    for pos in spots:
      x = pos[0]*pixel_sz - xbeam; y = pos[1]*pixel_sz - ybeam
      distance = sqrt(x*x+y*y)
      resolutions.append(
        info.C.lambda0*1E10 /(2.*sin(0.5*atan(distance/xtd)))
        )

    relevant_resolution = min(resolutions)

    recommended_sampling = relevant_resolution / largest_p1_cell_dimension

    recommended_sampling*=0.5  # 1/2 interval makes indexing much more reliable
    return recommended_sampling
