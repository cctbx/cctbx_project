from cctbx.array_family import flex
from cctbx.development import random_structure
from smtbx.refinement import constraints
from smtbx.refinement.constraints import adp as adp_constraints
from smtbx.refinement.constraints import geometrical_hydrogens

class random_xray_structure(random_structure.xray_structure):

  def __init__(self, space_group_info, u_iso_xor_u_aniso=True, **kwds):
    super(random_xray_structure, self).__init__(space_group_info, **kwds)
    if u_iso_xor_u_aniso: return
    if kwds['use_u_iso'] and kwds['use_u_aniso']:
      for sc in self.scatterers():
        sc.flags.set_use_u_iso(True).set_use_u_aniso(True)


class test_case(object):

  class __metaclass__(type):
    def __init__(cls, classname, bases, classdict):
      exercises = []
      for base in bases:
        try:
          exercises.extend(base.exercises)
        except AttributeError:
          pass
      for name, attr in classdict.items():
        if callable(attr) and name.startswith('exercise'):
          exercises.append(attr)
      dsu = [ (ex.__name__, ex) for ex in exercises ]
      dsu.sort()
      cls.exercises = [ ex for foo, ex in dsu ]

  def run(cls, verbose=False, *args, **kwds):
    if verbose: print cls.__name__
    for exercise in cls.exercises:
      if verbose: print "\t%s ... " % exercise.__name__,
      o = cls(*args, **kwds)
      exercise(o)
      if verbose: print "OK"
  run = classmethod(run)

def generate_hydrogen_constraints(structure, connectivity_table):
  # This is purely for testing purposes and is not in any way intended
  # to be a complete or comprehensive generation of constraints

  sc = structure.scatterers()
  conformer_indices = connectivity_table.conformer_indices
  sym_excl_indices = connectivity_table.sym_excl_indices
  sc_types = sc.extract_scattering_types()
  pair_sym_table = connectivity_table.pair_asu_table.extract_pair_sym_table(
    skip_j_seq_less_than_i_seq=False,
    all_interactions_from_inside_asu=True)
  h_constraints = []
  for i_seq, j_seq_dict in enumerate(pair_sym_table):
    conformer_i = conformer_indices[i_seq]
    if sc_types[i_seq] in ('H', 'D'): continue
    h_count = 0
    constrained_site_indices = []
    pivot_neighbour_count = 0
    for j_seq, sym_ops in j_seq_dict.items():
      if not (   conformer_i == 0
              or conformer_indices[j_seq] == 0
              or conformer_i == conformer_indices[j_seq]):
        continue
      if sc_types[j_seq] in ('H', 'D'):
        h_count += sym_ops.size()
        constrained_site_indices.append(j_seq)
      else: pivot_neighbour_count += sym_ops.size()
    rotating = False
    stretching = False
    constraint_type = None
    u_eq_multiplier = 1.2
    if h_count == 0: continue
    elif h_count == 1:
      if pivot_neighbour_count == 1:
        if sc_types[i_seq] == 'C':
          constraint_type = geometrical_hydrogens.terminal_linear_ch_site
        else:
          constraint_type = geometrical_hydrogens.terminal_tetrahedral_xh_site
          if sc_types[i_seq] == 'O':
            u_eq_multiplier  = 1.5
      elif pivot_neighbour_count == 2 and sc_types[i_seq] in ('C', 'N'):
        constraint_type = geometrical_hydrogens.secondary_planar_xh_site
        stretching = True
      elif pivot_neighbour_count == 3 and sc_types[i_seq] == 'C':
        constraint_type = geometrical_hydrogens.tertiary_ch_site
    elif h_count == 2:
      if pivot_neighbour_count == 1 and sc_types[i_seq] in ('C', 'N'):
        constraint_type = geometrical_hydrogens.terminal_planar_xh2_sites
      elif pivot_neighbour_count == 2 and sc_types[i_seq] == 'C':
        constraint_type = geometrical_hydrogens.secondary_ch2_sites
      elif pivot_neighbour_count == 0 and sc_types[i_seq] == 'O': # water
        u_eq_multiplier = 1.5
        for idx in constrained_site_indices:
          h_constraints.append(
            adp_constraints.u_iso_proportional_to_pivot_u_eq(
              u_eq_scatterer_idx=i_seq,
              u_iso_scatterer_idx=idx,
              multiplier=u_eq_multiplier))
    elif h_count == 3:
      if pivot_neighbour_count == 1:
        constraint_type = geometrical_hydrogens.terminal_tetrahedral_xh3_sites
        u_eq_multiplier = 1.5
    if constraint_type is not None:
      current = constraint_type(
        rotating=rotating,
        stretching=stretching,
        pivot=i_seq,
        constrained_site_indices=constrained_site_indices)
      h_constraints.append(current)
      if sc[i_seq].flags.use_u_iso(): continue # XXX u_iso can't yet ride on u_iso
      for idx in constrained_site_indices:
        h_constraints.append(
          adp_constraints.u_iso_proportional_to_pivot_u_eq(
            u_eq_scatterer_idx=i_seq,
            u_iso_scatterer_idx=idx,
            multiplier=u_eq_multiplier))
  return h_constraints
