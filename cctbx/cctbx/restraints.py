import cctbx.crystal.direct_space_asu
import cctbx.array_family.flex
import scitbx.array_family.shared

import boost.python
ext = boost.python.import_ext("cctbx_restraints_ext")
from cctbx_restraints_ext import *

import cctbx.sgtbx
import scitbx.stl.map
import scitbx.stl.set
import scitbx.stl.vector

pair_sym_ops = cctbx.sgtbx.stl_vector_rt_mx

pair_asu_j_sym_groups = scitbx.stl.vector.set_unsigned
pair_asu_j_sym_group = scitbx.stl.set.unsigned

repulsion_radius_table = scitbx.stl.map.stl_string_double

repulsion_distance_table = scitbx.stl.map.stl_string_stl_map_stl_string_double
repulsion_distance_dict = scitbx.stl.map.stl_string_double
