# LIBTBX_SET_DISPATCHER_NAME cxi.stream_to_pickle
""" Utility for converting stream files from CrystFEL version 0.5.3 to
cctbx.xfel pickle files.
"""

from __future__ import absolute_import, division, print_function
import re
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell
from cctbx.crystal import symmetry
from cctbx import miller
import logging

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

EV_PER_A = 12398.4187
# Regular expressions, set up so one can use groups to extract the data
re_energy = re.compile("photon_energy_eV\s=\s([0-9]+\.[0-9]+)")
re_uc = re.compile("""Cell\sparameters\s
                   ([0-9]+\.[0-9]+)\s # a
                   ([0-9]+\.[0-9]+)\s # b
                   ([0-9]+\.[0-9]+)\snm,\s # c
                   ([0-9]+\.[0-9]+)\s # alpha
                   ([0-9]+\.[0-9]+)\s # beta
                   ([0-9]+\.[0-9]+) # gamma""", re.X)

# groups 1-7: h,k,l,I,sig(i),peak, background,fs,ss
re_miller = re.compile("""\s*(-?[0-9]{1,3})
                          \s*(-?[0-9]{1,3})
                          \s*(-?[0-9]{1,3}) #h,k,l
                          \s*(-?[0-9]+\.[0-9]+)
                          \s*(-?[0-9]+\.[0-9]+) #I, SigI
                          ([\s*-?[0-9]+\.[0-9]+){2}  #Peak, background
                          \s*([0-9]+\.[0-9]+)
                          \s*([0-9]+\.[0-9]+) #fs, ssi""", re.X)
re_lattice_type = re.compile("lattice_type\s=\s([a-zA-Z]+)")
re_centering = re.compile("centering\s=\s([A-Z])")

# note the setting lattice is in nm^-1
re_Astar = re.compile("""astar\s=\s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)""", re.X)
re_Bstar = re.compile("""bstar\s=\s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)""", re.X)
re_Cstar = re.compile("""cstar\s=\s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)
                         \s*(-?\+?[0-9]+\.[0-9]+)""", re.X)


def image_template():
    return {'Millers': [], 'Is': [], 'sigIs': [], 'location': []}


def check_image(image):
  if things_in_image == set(image):
    if len(image['Millers']) > 2:
      return True
  else:
    return False


def crystfel_to_cctbx_coord_system(abasis, bbasis, cbasis):
  """CrystFEL is RHS with z down the beam, and y to the ceiling.
  CCTBX.Postrefine is RHS with z to the source, and y to the ceiling.
  """
  from scitbx.matrix import sqr
  from cctbx import crystal_orientation

  a_mat = sqr((abasis[0], bbasis[0], cbasis[0],
               abasis[1], bbasis[1], cbasis[1],
               abasis[2], bbasis[2], cbasis[2]))

  coord_transformation = sqr((  -1,  0,   0,
                                 0,  1,   0,
                                 0,  0,  -1))
  new_coords = a_mat.__mul__(coord_transformation)

  ori = crystal_orientation.crystal_orientation(new_coords, crystal_orientation.basis_type.reciprocal)
  logging.debug("\naStar: {}\nbStar: {}\ncStar: {}".format(abasis, bbasis, cbasis))
  logging.debug(str(ori))
  return ori


def unit_cell_to_symetry_object(img_dict):
  xsym = symmetry(unit_cell=img_dict['current_orientation'][0].unit_cell().niggli_cell(),
                  space_group_symbol=point_group)
  miller_set = miller.set(crystal_symmetry=xsym,
                          indices=flex.miller_index(img_dict['Millers']),
                          anomalous_flag=True)
  miller_array = miller.array(miller_set,
                              flex.double(img_dict['Is']),
                              flex.double(img_dict['sigIs']))
  miller_array.set_observation_type_xray_intensity()
  miller_array.set_info("Raw data obtained by integration using CrystFEL")
  return miller_array


def make_int_pickle(img_dict, filename):
  img_dict['current_orientation'] = [crystfel_to_cctbx_coord_system(
                    img_dict['aStar'],
                    img_dict['bStar'],
                    img_dict['cStar'],)]
  try:
    final_dict = {'observations': [unit_cell_to_symetry_object(img_dict)],
                  'mapped_predictions': [flex.vec2_double(img_dict['location'])],
                  'xbeam': 96.99,
                  'ybeam': 96.97,
                  'distance': 150.9,
                  "sa_parameters": ['None'],
                  "pointgroup": point_group,
                  "unit_cell": img_dict['unit cell'],
                  'current_orientation': img_dict['current_orientation'],
                  'wavelength': img_dict['wavelength'],
                  'centering': img_dict['centering'],
                  'lattice_type': img_dict['lattice_type']}
    pickle.dump(final_dict, open(filename, 'wb'))
    logging.info("dumped image {}".format(filename))
  except AssertionError as a:
    logging.warning("Failed an assertion on image {}! {}".format(filename, a.message))
  except Exception:
    logging.warning("Failed to make a dictionairy for image {}".format(filename))

if __name__ == "__main__":

  logging.critical("NOT READY FOR PRIME-TIME. CONTACT ZELDIN@STANFORD.EDU IF YOU WANT TO USE THIS.")

  import argparse
  parser = argparse.ArgumentParser(description=
                                 ('Create indexing pickles from a'
                                  'crystfel stream file.'))
  parser.add_argument('filename', type=str, nargs=1,
                    help='The filename of the stream file to be converted.')
  parser.add_argument('point_group', type=str, nargs=1,
                    help='The space group to be assigned.')
  parser.add_argument('--tag', type=str, nargs=1,
                    help="Prefix to be used for the indexing pickles.")
  args = parser.parse_args()
  stream_file = args.filename[0]
  point_group = args.point_group[0]
  if args.tag:
    tag = args.tag[0]
  else:
    tag = stream_file[0].split('.')[0]


  things_in_image = {'Millers', 'Is', 'sigIs', 'unit cell', 'aStar', 'bStar',
                     'cStar', 'wavelength', 'centering', 'lattice_type',
                     'location'}
  with open(stream_file, "r+") as stream:
    count = 1
    this_image = image_template()
    for line in stream:
      millers = re_miller.match(line)
      if millers:
        this_image["Millers"].append((int(millers.group(1)),
                                      int(millers.group(2)),
                                      int(millers.group(3))))
        this_image["Is"].append(float(millers.group(4)))
        this_image["sigIs"].append(float(millers.group(5)))
        this_image["location"].append((float(millers.group(7)),
                                       float(millers.group(8))))
        continue

      uc = re_uc.match(line)
      if uc:  # Start of a new crystal dictionary
        if check_image(this_image):  # i.e. it's a complete image dictionairy
          make_int_pickle(this_image, "{}_{:04d}.pickle".format(tag, count))
          count += 1
          this_image = image_template()
          this_image["unit cell"] = unit_cell((float(uc.group(1)) * 10,  # nm to A
                                               float(uc.group(2)) * 10,
                                               float(uc.group(3)) * 10,
                                               float(uc.group(4)),
                                               float(uc.group(5)),
                                               float(uc.group(6))))

        else:
          this_image = image_template()
          this_image["unit cell"] = unit_cell((float(uc.group(1)) * 10,  # nm to A
                                               float(uc.group(2)) * 10,
                                               float(uc.group(3)) * 10,
                                               float(uc.group(4)),
                                               float(uc.group(5)),
                                               float(uc.group(6))))
        continue

      energy = re_energy.match(line)
      if energy:
        this_image['wavelength'] = float(energy.group(1)) / EV_PER_A
        continue

      lattice = re_lattice_type.match(line)
      if lattice:
        this_image['lattice_type'] = lattice.group(1)
        continue

      centering = re_centering.match(line)
      if centering:
        this_image['centering'] = centering.group(1)
        continue

      astar = re_Astar.match(line)
      if astar:
        this_image['aStar'] = [float(astar.group(1)) / 10,
                               float(astar.group(2)) / 10,
                               float(astar.group(3)) / 10]
        continue

      bstar = re_Bstar.match(line)
      if bstar:
        this_image['bStar'] = [float(bstar.group(1)) / 10,
                               float(bstar.group(2)) / 10,
                               float(bstar.group(3)) / 10]
        continue

      cstar = re_Cstar.match(line)
      if cstar:
        this_image['cStar'] = [float(cstar.group(1)) / 10,
                               float(cstar.group(2)) / 10,
                               float(cstar.group(3)) / 10]

        continue
# After this for loop, we should have an array of dictionairies,
# each containing the info needed.
# Just need to pickle the array of dictionairies! et voila!
