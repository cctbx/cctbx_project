from __future__ import absolute_import, division, print_function
from StringIO import StringIO
# import sys, time

# from libtbx.utils import Sorry
# from scitbx.array_family import flex
# from libtbx import adopt_init_args
# from libtbx.str_utils import make_header

from cctbx.geometry_restraints.manager import manager as standard_manager

orca_master_phil_str = '''
  use_orca = True
    .type = bool
  method = PM6
    .type = str
  basis_set = None
    .type = str
  solvent_model = CPCM
    .type = str
'''

class manager(standard_manager):
  def __init__(self,
               # pdb_hierarchy,
               params,
               energy_components=None,
               # gradients=None,
               # gradients_factory=None, #flex.vec3_double,
               log=StringIO()):
    # self.gradients_factory = gradients_factory
    adopt_init_args(self, locals(), exclude=["log"])

def digester(standard_geometry_restraints_manager,
             params,
             log=StringIO(),
             ):
  sgrm = standard_geometry_restraints_manager
  qi_grm = manager(params, log=log)
  for attr, value in vars(sgrm).items(): setattr(qi_grm, attr, value)
  qi_grm.standard_geometry_restraints_manager = sgrm
  return qi_grm
