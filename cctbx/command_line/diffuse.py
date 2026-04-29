"""simulate diffuse scattering"""
# LIBTBX_SET_DISPATCHER_NAME phenix.diffuse

from __future__ import absolute_import, division, print_function
import sys
import iotbx.pdb
import math
from cctbx.array_family import flex
from cctbx import miller
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from six.moves import range

msg="""
Description
-----------
This program calculates reciprocal space diffuse scattering maps from
multi-model PDB ensembles.

Output
------
MTZ file containing diffuse scattering intensity values (I) at each Bragg peak.

Command line Usage
------------------
phenix.diffuse pdb=test.pdb probabilities=0.3,0.3,0.3,0.1 resolution=4.0 prefix=tst

pdb: input PDB file. CRYST1 symmetry description is required.

probabilities: the weighted probability for each model in the corresponding PDB
file. The values listed in this option will be assigned to the PDB models
starting with MODEL 0. Alternatively, leaving out the probabilities option will
lead to equal weights for each model in the PDB file.

resolution: the desired d_min for the output MTZ data set

prefix: the desired MTZ file name

Author
------
Andrew Van Benschoten (andrew.vanbenschoten@ucsf.edu)
"""

def run(arg):
  args = get_input_dict(arg)
  if(len(args)!=4):
    msg="""Bad inputs.
Usage example:
  phenix.diffuse pdb=m.pdb probabilities=0.5,0.5 resolution=4.0 prefix=tst"""
    raise Sorry(msg)
  data = ensemble(
    pdb_file_name = args['pdb'],
    probabilities = args['probabilities'])
  data.get_models()
  for model in data.models:
    model.initialize(d_min = float(args['resolution']))
  diffuse(
    models           = data.models,
    crystal_symmetry = data.symmetry,
    scale_factor     = 1).write_mtz_file(prefix = args['prefix'])

class ensemble(object):

  def __init__(self, pdb_file_name, probabilities):
    adopt_init_args(self, locals())
    pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    self.hierarchy = pdb_inp.construct_hierarchy()
    self.symmetry = pdb_inp.crystal_symmetry_from_cryst1()
    self.xray_structures_p1 = [xrs.expand_to_p1() for xrs in
      pdb_inp.xray_structures_simple(crystal_symmetry=self.symmetry)]

  def get_models(self):
    self.models = []
    models = self.hierarchy.models()
    weights = []
    if(self.probabilities is not None):
      new_probs = get_probabilities(self.probabilities)
      for i in range(0,len(new_probs)):
        d = float(new_probs[i])
        weights.append(d)
      if(len(new_probs) != len(models)):
        raise Sorry(
          "The number of models and number of given probabilities must match")
    else:
      for model_ in models:
        d = float(1/len(models))
        weights.append(d)
    for i, model_ in enumerate(models):
      m = model(
        model             = model,
        xray_structure_p1 = self.xray_structures_p1[i],
        probability       = weights[i])
      self.models.append(m)

class model(object):

  def __init__(self, model, xray_structure_p1, probability):
    adopt_init_args(self, locals())
    assert xray_structure_p1.crystal_symmetry().space_group().type().number()==1

  def initialize(self, d_min):
    f = self.xray_structure_p1.structure_factors(d_min=d_min).f_calc()
    self.f_weighted = f*self.probability
    f_squared = abs(f).set_observation_type_xray_amplitude().f_as_f_sq()
    self.f_squared_weighted = f_squared*self.probability

class diffuse(object):
  "Class for all diffuse maps produced in reciprocal space"

  #REMOVE SAMPLING
  def __init__(self, models, crystal_symmetry, scale_factor):
    adopt_init_args(self, locals())
    self.lattice = {}
    self.calculate_map()

  def calculate_map(self):
    sum_fc = None
    sum_fc_square = None
    for model in self.models:
      x, y = model.f_weighted.data(), model.f_squared_weighted.data()
      if sum_fc is None:
        sum_fc = x
        sum_fc_square = y
      else:
        sum_fc = sum_fc + x
        sum_fc_square = sum_fc_square + y
    self.diffuse_signal = self.models[0].f_weighted.customized_copy(
      data = sum_fc_square - flex.abs(sum_fc)**2)
    self.write_squared_amplitudes(miller_array = self.diffuse_signal)

  def write_squared_amplitudes(self, miller_array):
    eps = 1.e-9
    for hkl, intensity in miller_array:
      h_int = hkl[0]
      k_int = hkl[1]
      l_int = hkl[2]
      intensity_new = intensity/self.scale_factor
      if h_int not in self.lattice:
        self.lattice[h_int] = {}
      if k_int not in self.lattice[h_int]:
        self.lattice[h_int][k_int] = {}
      if l_int not in self.lattice[h_int][k_int]:
        self.lattice[h_int][k_int][l_int] = 0
      self.lattice[h_int][k_int][l_int] += intensity_new

  def write_mtz_file(self, prefix):
    indices = flex.miller_index()
    i_obs = flex.double()
    sig_i = flex.double()
    assert self.scale_factor != 0
    for key_h in self.lattice:
      for key_k in self.lattice[key_h]:
        for key_l in self.lattice[key_h][key_k]:
          indices.append([key_h, key_k, key_l])
          # why convert to int and then go back to float? must be a bug?..
          io = float("%4d"%self.lattice[key_h][key_k][key_l])/self.scale_factor
          i_obs.append(io)
          sig_i.append(math.sqrt(io))
    # get miller array object
    ma = miller.array(miller_set=miller.set(self.crystal_symmetry, indices),
      data=i_obs, sigmas=sig_i)
    ma.set_observation_type_xray_intensity()
    mtz_dataset = ma.as_mtz_dataset(column_root_label="I")
    mtz_dataset.mtz_object().write(prefix + '.mtz')

def get_input_dict(args):
  dic = dict()
  for arg in args:
    spl=arg.split('=')
    if len(spl)==2:
      dic[spl[0]] = spl[1]
  if 'probabilities' not in dic:
    dic['probabilities'] = None
  return dic

def get_probabilities(input):
  result = flex.double([float(d) for d in input.split(',')])
  if(abs(1.0-flex.sum(result))>1.e-3):
    raise Sorry("Sorry, the given probabilities must sum to one")
  return result

if __name__ == '__main__':
  run(sys.argv[1:])

