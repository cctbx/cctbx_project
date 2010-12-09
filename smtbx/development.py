from cctbx.development import random_structure
from cctbx import crystal, xray
from cctbx.array_family import flex
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

def sucrose():
  return xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(7.783, 8.7364, 10.9002, 90, 102.984, 90),
      space_group_symbol='hall:  P 2yb'),
    scatterers=flex.xray_scatterer((
      xray.scatterer( #0
                      label='O1',
                      site=(-0.131694, 0.935444, 0.877178),
                      u=(0.000411, 0.000231, 0.000177,
                         -0.000017, 0.000112, -0.000003)),
      xray.scatterer( #1
                      label='O2',
                      site=(-0.214093, 0.788106, 1.081133),
                      u=(0.000809, 0.000424, 0.000219,
                       0.000105, 0.000177, 0.000051)),
      xray.scatterer( #2
                      label='H2',
                      site=(-0.227111, 0.876745, 1.102068),
                      u=0.050326),
      xray.scatterer( #3
                      label='O3',
                      site=(-0.145099, 0.520224, 0.848374),
                      u=(0.000860, 0.000240, 0.000587,
                         -0.000149, 0.000303, -0.000093)),
      xray.scatterer( #4
                      label='H3',
                      site=(-0.073357, 0.454114, 0.840907),
                      u=0.064297),
      xray.scatterer( #5
                      label='O4',
                      site=(0.202446, 0.586695, 0.808946),
                      u=(0.000848, 0.000352, 0.000315,
                         0.000294, 0.000256, 0.000111)),
      xray.scatterer( #6
                      label='H4',
                      site=(0.234537, 0.569702, 0.743559),
                      u=0.052973),
      xray.scatterer( #7
                      label='O5',
                      site=(0.247688, 0.897720, 0.729153),
                      u=(0.000379, 0.000386, 0.000270,
                         -0.000019, 0.000106, 0.000001)),
      xray.scatterer( #8
                      label='H5',
                      site=(0.323942, 0.957695, 0.764491),
                      u=0.040241),
      xray.scatterer( #9
                      label='O6',
                      site=(-0.183682, 1.239489, 0.712253),
                      u=(0.000240, 0.000270, 0.000208,
                         0.000035, 0.000002, -0.000057)),
      xray.scatterer( #10
                      label='O7',
                      site=(-0.460944, 1.095427, 0.826725),
                      u=(0.000454, 0.000438, 0.000330,
                         0.000029, 0.000103, 0.000091)),
      xray.scatterer( #11
                      label='H7',
                      site=(-0.360598, 1.061356, 0.849066),
                      u=0.048134),
      xray.scatterer( #12
                      label='O8',
                      site=(-0.590098, 1.234857, 0.477504),
                      u=(0.000341, 0.000376, 0.000283,
                         0.000054, -0.000022, 0.000051)),
      xray.scatterer( #13
                      label='H8',
                      site=(-0.596179, 1.326096, 0.493669),
                      u=0.041878),
      xray.scatterer( #14
                      label='O9',
                      site=(-0.294970, 1.016028, 0.425505),
                      u=(0.000509, 0.000266, 0.000191,
                         -0.000031, 0.000056, -0.000050)),
      xray.scatterer( #15
                      label='H9',
                      site=(-0.305666, 0.937198, 0.463939),
                      u=0.035859),
      xray.scatterer( #16
                      label='O10',
                      site=(0.121231, 1.098881, 0.529558),
                      u=(0.000428, 0.000442, 0.000264,
                         0.000037, 0.000149, 0.000061)),
      xray.scatterer( #17
                      label='H10',
                      site=(0.164296, 1.026179, 0.573427),
                      u=0.042667),
      xray.scatterer( #18
                      label='O11',
                      site=(-0.108540, 0.986958, 0.671415),
                      u=(0.000301, 0.000160, 0.000170,
                         -0.000018, 0.000028, 0.000003)),
      xray.scatterer( #19
                      label='C1',
                      site=(-0.205698, 0.782086, 0.859435),
                      u=(0.000378, 0.000249, 0.000211,
                         -0.000046, 0.000020, 0.000007)),
      xray.scatterer( #20
                      label='H1',
                      site=(-0.282033, 0.773748, 0.774955),
                      u=0.033128),
      xray.scatterer( #21
                      label='C2',
                      site=(-0.315620, 0.763019, 0.956870),
                      u=(0.000478, 0.000387, 0.000268,
                         0.000006, 0.000115, 0.000085)),
      xray.scatterer( #22
                      label='H2B',
                      site=(-0.413052, 0.834909, 0.939327),
                      u=0.042991),
      xray.scatterer( #23
                      label='C3',
                      site=(-0.057676, 0.663433, 0.874279),
                      u=(0.000481, 0.000203, 0.000211,
                         -0.000039, 0.000097, -0.000026)),
      xray.scatterer( #24
                      label='H3B',
                      site=(0.010861, 0.664198, 0.961490),
                      u=0.033010),
      xray.scatterer( #25
                      label='C4',
                      site=(0.064486, 0.697929, 0.785789),
                      u=(0.000516, 0.000251, 0.000171,
                         0.000108, 0.000066, 0.000006)),
      xray.scatterer( #26
                      label='H2A',
                      site=(-0.364431, 0.660423, 0.951038),
                      u=0.042991),
      xray.scatterer( #27
                      label='H4B',
                    site=(-0.001643, 0.690588, 0.698163),
                    u=0.034151),
      xray.scatterer( #28
                      label='C5',
                      site=(0.134552, 0.859023, 0.812706),
                      u=(0.000346, 0.000319, 0.000114,
                         0.000016, 0.000027, 0.000034)),
      xray.scatterer( #29
                      label='H5B',
                      site=(0.204847, 0.862433, 0.899373),
                      u=0.028874),
      xray.scatterer( #30
                      label='C6',
                      site=(-0.013996, 0.975525, 0.800238),
                      u=(0.000321, 0.000194, 0.000125,
                         -0.000031, 0.000037, -0.000006)),
      xray.scatterer( #31
                      label='H6',
                      site=(0.037624, 1.075652, 0.827312),
                      u=0.023850),
      xray.scatterer( #32
                      label='C7',
                      site=(-0.130206, 1.141642, 0.624140),
                      u=(0.000313, 0.000143, 0.000184,
                         -0.000000, 0.000045, -0.000001)),
      xray.scatterer( #33
                      label='C8',
                      site=(0.043311, 1.203061, 0.603290),
                      u=(0.000354, 0.000269, 0.000247,
                         -0.000020, 0.000085, 0.000024)),
      xray.scatterer( #34
                      label='H8A',
                      site=(0.023319, 1.300964, 0.560480),
                      u=0.034035),
      xray.scatterer( #35
                      label='H8B',
                      site=(0.123967, 1.219261, 0.684061),
                      u=0.034035),
      xray.scatterer( #36
                      label='C9',
                      site=(-0.285452, 1.143037, 0.507357),
                      u=(0.000294, 0.000203, 0.000181,
                         0.000009, 0.000062, 0.000042)),
      xray.scatterer( #37
                      label='H9B',
                      site=(-0.273079, 1.235040, 0.458631),
                      u=0.026214),
      xray.scatterer( #38
                      label='C10',
                      site=(-0.444628, 1.167026, 0.565148),
                      u=(0.000289, 0.000234, 0.000210,
                         0.000018, 0.000036, 0.000040)),
      xray.scatterer( #39
                      label='H10B',
                      site=(-0.480820, 1.069195, 0.595409),
                      u=0.029493),
      xray.scatterer( #40
                      label='C11',
                      site=(-0.371995, 1.272625, 0.676806),
                      u=(0.000324, 0.000199, 0.000264,
                         0.000054, 0.000089, 0.000009)),
      xray.scatterer( #41
                      label='H11',
                      site=(-0.388144, 1.379127, 0.648319),
                      u=0.031434),
      xray.scatterer( #42
                      label='C12',
                      site=(-0.453277, 1.252159, 0.788677),
                      u=(0.000396, 0.000422, 0.000263,
                         0.000043, 0.000125, -0.000061)),
      xray.scatterer( #43
                      label='H12B',
                      site=(-0.571879, 1.293694, 0.768457),
                      u=0.041343),
      xray.scatterer( #44
                      label='H12A',
                      site=(-0.385488, 1.310416, 0.858834),
                      u=0.041343)
    )))
