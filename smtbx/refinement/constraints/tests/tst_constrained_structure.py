from cctbx import crystal, xray
from cctbx.array_family import flex
from smtbx.refinement import constraints
import smtbx.refinement.constraints.geometrical_hydrogens as _
import smtbx.refinement.constraints as core
import itertools

def exercise():
  # sucrose from Olex 2 samples
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(7.783, 8.7364, 10.9002, 90, 102.984, 90),
      space_group_symbol='hall:  P 2yb'),
    scatterers=flex.xray_scatterer((
      xray.scatterer(
        label='O1',
        site=(-0.131694, 0.935444, 0.877178),
        u=(0.000411, 0.000231, 0.000177,
           -0.000017, 0.000112, -0.000003)),
      xray.scatterer(
        label='O2',
        site=(-0.214093, 0.788106, 1.081133),
        u=(0.000809, 0.000424, 0.000219,
           0.000105, 0.000177, 0.000051)),
      xray.scatterer(
        label='H2',
        site=(-0.227111, 0.876745, 1.102068),
        u=0.050326),
      xray.scatterer(
        label='O3',
        site=(-0.145099, 0.520224, 0.848374),
        u=(0.000860, 0.000240, 0.000587,
           -0.000149, 0.000303, -0.000093)),
      xray.scatterer(
        label='H3'),
      xray.scatterer(
        label='O4',
        site=(0.202446, 0.586695, 0.808946),
        u=(0.000848, 0.000352, 0.000315,
           0.000294, 0.000256, 0.000111)),
      xray.scatterer(
        label='H4'),
      xray.scatterer(
        label='O5',
        site=(0.247688, 0.897720, 0.729153),
        u=(0.000379, 0.000386, 0.000270,
           -0.000019, 0.000106, 0.000001)),
      xray.scatterer(
        label='H5'),
      xray.scatterer(
        label='O6',
        site=(-0.183682, 1.239489, 0.712253),
        u=(0.000240, 0.000270, 0.000208,
           0.000035, 0.000002, -0.000057)),
      xray.scatterer(
        label='O7',
        site=(-0.460944, 1.095427, 0.826725),
        u=(0.000454, 0.000438, 0.000330,
           0.000029, 0.000103, 0.000091)),
      xray.scatterer(
        label='H7'),
      xray.scatterer(
        label='O8',
        site=(-0.590098, 1.234857, 0.477504),
        u=(0.000341, 0.000376, 0.000283,
           0.000054, -0.000022, 0.000051)),
      xray.scatterer(
        label='H8'),
      xray.scatterer(
        label='O9',
        site=(-0.294970, 1.016028, 0.425505),
        u=(0.000509, 0.000266, 0.000191,
           -0.000031, 0.000056, -0.000050)),
      xray.scatterer(
        label='H9'),
      xray.scatterer(
        label='O10',
        site=(0.121231, 1.098881, 0.529558),
        u=(0.000428, 0.000442, 0.000264,
           0.000037, 0.000149, 0.000061)),
      xray.scatterer(
        label='H10'),
      xray.scatterer(
        label='O11',
        site=(-0.108540, 0.986958, 0.671415),
        u=(0.000301, 0.000160, 0.000170,
           -0.000018, 0.000028, 0.000003)),
      xray.scatterer(
        label='C1',
        site=(-0.205698, 0.782086, 0.859435),
        u=(0.000378, 0.000249, 0.000211,
           -0.000046, 0.000020, 0.000007)),
      xray.scatterer(
        label='H1'),
      xray.scatterer(
        label='C2',
        site=(-0.315620, 0.763019, 0.956870),
        u=(0.000478, 0.000387, 0.000268,
           0.000006, 0.000115, 0.000085)),
      xray.scatterer(
        label='H2B',
        site=(-0.413052, 0.834909, 0.939327),
        u=0.042991),
      xray.scatterer(
        label='H2A'),
      xray.scatterer(
        label='C3',
        site=(-0.057676, 0.663433, 0.874279),
        u=(0.000481, 0.000203, 0.000211,
           -0.000039, 0.000097, -0.000026)),
      xray.scatterer(
        label='H3B'),
      xray.scatterer(
        label='C4',
        site=(0.064486, 0.697929, 0.785789),
        u=(0.000516, 0.000251, 0.000171,
           0.000108, 0.000066, 0.000006)),
      xray.scatterer(
        label='H4B'),
      xray.scatterer(
        label='C5',
        site=(0.134552, 0.859023, 0.812706),
        u=(0.000346, 0.000319, 0.000114,
           0.000016, 0.000027, 0.000034)),
      xray.scatterer(
        label='H5B'),
      xray.scatterer(
        label='C6',
        site=(-0.013996, 0.975525, 0.800238),
        u=(0.000321, 0.000194, 0.000125,
           -0.000031, 0.000037, -0.000006)),
      xray.scatterer(
        label='H6'),
      xray.scatterer(
        label='C7',
        site=(-0.130206, 1.141642, 0.624140),
        u=(0.000313, 0.000143, 0.000184,
           -0.000000, 0.000045, -0.000001)),
      xray.scatterer(
        label='C8',
        site=(0.043311, 1.203061, 0.603290),
        u=(0.000354, 0.000269, 0.000247,
           -0.000020, 0.000085, 0.000024)),
      xray.scatterer(
        label='H8A'),
      xray.scatterer(
        label='H8B'),
      xray.scatterer(
        label='C9',
        site=(-0.285452, 1.143037, 0.507357),
        u=(0.000294, 0.000203, 0.000181,
           0.000009, 0.000062, 0.000042)),
      xray.scatterer(
        label='H9B'),
      xray.scatterer(
        label='C10',
        site=(-0.444628, 1.167026, 0.565148),
        u=(0.000289, 0.000234, 0.000210,
           0.000018, 0.000036, 0.000040)),
      xray.scatterer(
        label='H10B'),
      xray.scatterer(
        label='C11',
        site=(-0.371995, 1.272625, 0.676806),
        u=(0.000324, 0.000199, 0.000264,
           0.000054, 0.000089, 0.000009)),
      xray.scatterer(
        label='H11'),
      xray.scatterer(
        label='C12',
        site=(-0.453277, 1.252159, 0.788677),
        u=(0.000396, 0.000422, 0.000263,
           0.000043, 0.000125, -0.000061)),
      xray.scatterer(
        label='H12B'),
      xray.scatterer(
        label='H12A')
    )))
  for sc in xs.scatterers():
    sc.flags.set_grad_site(True)
    sc.flags.set_grad_u_aniso(True)

  reparametrisation = constraints.reparametrisation(
    structure=xs,
    geometrical_constraints=[
      _.terminal_tetrahedral_xh_site(
        pivot=1,
        constrained_site_indices=(2,)),
      _.terminal_tetrahedral_xh_site(
        pivot=3,
        constrained_site_indices=(4,)),
      _.terminal_tetrahedral_xh_site(
        pivot=5,
        constrained_site_indices=(6,)),
      _.terminal_tetrahedral_xh_site(
        pivot=7,
        constrained_site_indices=(8,)),
      _.terminal_tetrahedral_xh_site(
        pivot=10,
        constrained_site_indices=(11,)),
      _.terminal_tetrahedral_xh_site(
        pivot=12,
        constrained_site_indices=(13,)),
      _.terminal_tetrahedral_xh_site(
        pivot=14,
        constrained_site_indices=(15,)),
      _.terminal_tetrahedral_xh_site(
        pivot=16,
        constrained_site_indices=(17,)),
      _.tertiary_ch_site(
        pivot=19,
        constrained_site_indices=(20,)),
      _.secondary_ch2_sites(
        pivot=21,
        constrained_site_indices=(22, 23)),
      _.tertiary_ch_site(
        pivot=24,
        constrained_site_indices=(25,)),
      _.tertiary_ch_site(
        pivot=26,
        constrained_site_indices=(27,)),
      _.tertiary_ch_site(
        pivot=28,
        constrained_site_indices=(29,)),
      _.tertiary_ch_site(
        pivot=30,
        constrained_site_indices=(31,)),
      _.secondary_ch2_sites(
        pivot=33,
        constrained_site_indices=(34, 35)),
      _.tertiary_ch_site(
        pivot=36,
        constrained_site_indices=(37,)),
      _.tertiary_ch_site(
        pivot=38,
        constrained_site_indices=(39,)),
      _.tertiary_ch_site(
        pivot=40,
        constrained_site_indices=(41,)),
      _.secondary_ch2_sites(
        pivot=42,
        constrained_site_indices=(43, 44)),
      ])

  warned_once = False
  for sc, params in itertools.izip(reparametrisation.structure.scatterers(),
                                   reparametrisation.asu_scatterer_parameters):
    if sc.scattering_type != 'H':
      assert isinstance(params.site, core.independent_site_parameter)
      assert isinstance(params.u, core.independent_cartesian_adp)
      assert isinstance(params.occupancy, core.independent_scalar_parameter)
    else:
      expected = {
        "H2": (core.terminal_tetrahedral_xh_site, 'O2'),
        "H3": (core.terminal_tetrahedral_xh_site, 'O3'),
        "H4": (core.terminal_tetrahedral_xh_site, 'O4'),
        "H5": (core.terminal_tetrahedral_xh_site, 'O5'),
        "H7": (core.terminal_tetrahedral_xh_site, 'O7'),
        "H8": (core.terminal_tetrahedral_xh_site, 'O8'),
        "H9": (core.terminal_tetrahedral_xh_site, 'O9'),
        "H10": (core.terminal_tetrahedral_xh_site, 'O10'),
        "H1": (core.tertiary_ch_site, 'C1'),
        "H2A": (core.secondary_ch2_sites, 'C2'),
        "H2B": (core.secondary_ch2_sites, 'C2'),
        "H3B": (core.tertiary_ch_site, 'C3'),
        "H4B": (core.tertiary_ch_site, 'C4'),
        "H5B": (core.tertiary_ch_site, 'C5'),
        "H6": (core.tertiary_ch_site, 'C6'),
        "H8A": (core.secondary_ch2_sites, 'C8'),
        "H8B": (core.secondary_ch2_sites, 'C8'),
        "H9B": (core.tertiary_ch_site, 'C9'),
        "H10B": (core.tertiary_ch_site, 'C10'),
        "H11": (core.tertiary_ch_site, 'C11'),
        "H12A": (core.secondary_ch2_sites, 'C12'),
        "H12B": (core.secondary_ch2_sites, 'C12'),
        }.get(sc.label)
      if expected is None:
        if not warned_once:
          print "Warning: incomplete test coverage for H constraint types"
          warned_once = True
        continue
      expected_type, expected_pivot = expected
      assert isinstance(params.site, expected_type), \
             (sc.label, params.site, expected_type)
      assert ([ sc1.label for sc1 in params.site.argument(0).scatterers() ]
              == [expected_pivot]), sc.label

def run():
  import libtbx.utils
  libtbx.utils.show_times_at_exit()
  exercise()

if __name__ == '__main__':
  run()
