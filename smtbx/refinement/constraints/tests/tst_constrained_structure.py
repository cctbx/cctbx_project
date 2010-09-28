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
        label='H3'),
      xray.scatterer( #5
        label='O4',
        site=(0.202446, 0.586695, 0.808946),
        u=(0.000848, 0.000352, 0.000315,
           0.000294, 0.000256, 0.000111)),
      xray.scatterer( #6
        label='H4'),
      xray.scatterer( #7
        label='O5',
        site=(0.247688, 0.897720, 0.729153),
        u=(0.000379, 0.000386, 0.000270,
           -0.000019, 0.000106, 0.000001)),
      xray.scatterer( #8
        label='H5'),
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
        label='H7'),
      xray.scatterer( #12
        label='O8',
        site=(-0.590098, 1.234857, 0.477504),
        u=(0.000341, 0.000376, 0.000283,
           0.000054, -0.000022, 0.000051)),
      xray.scatterer( #13
        label='H8'),
      xray.scatterer( #14
        label='O9',
        site=(-0.294970, 1.016028, 0.425505),
        u=(0.000509, 0.000266, 0.000191,
           -0.000031, 0.000056, -0.000050)),
      xray.scatterer( #15
        label='H9'),
      xray.scatterer( #16
        label='O10',
        site=(0.121231, 1.098881, 0.529558),
        u=(0.000428, 0.000442, 0.000264,
           0.000037, 0.000149, 0.000061)),
      xray.scatterer( #17
        label='H10'),
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
        label='H1'),
      xray.scatterer( #21
        label='C2',
        site=(-0.315620, 0.763019, 0.956870),
        u=(0.000478, 0.000387, 0.000268,
           0.000006, 0.000115, 0.000085)),
      xray.scatterer( #22
        label='H2B'),
      xray.scatterer( #23
        label='C3',
        site=(-0.057676, 0.663433, 0.874279),
        u=(0.000481, 0.000203, 0.000211,
           -0.000039, 0.000097, -0.000026)),
      xray.scatterer( #24
        label='H3B'),
      xray.scatterer( #25
        label='C4',
        site=(0.064486, 0.697929, 0.785789),
        u=(0.000516, 0.000251, 0.000171,
           0.000108, 0.000066, 0.000006)),
      xray.scatterer( #26
        label='H2A'),
      xray.scatterer( #27
        label='H4B'),
      xray.scatterer( #28
        label='C5',
        site=(0.134552, 0.859023, 0.812706),
        u=(0.000346, 0.000319, 0.000114,
           0.000016, 0.000027, 0.000034)),
      xray.scatterer( #29
        label='H5B'),
      xray.scatterer( #30
        label='C6',
        site=(-0.013996, 0.975525, 0.800238),
        u=(0.000321, 0.000194, 0.000125,
           -0.000031, 0.000037, -0.000006)),
      xray.scatterer( #31
        label='H6'),
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
        label='H8A'),
      xray.scatterer( #35
        label='H8B'),
      xray.scatterer( #36
        label='C9',
        site=(-0.285452, 1.143037, 0.507357),
        u=(0.000294, 0.000203, 0.000181,
           0.000009, 0.000062, 0.000042)),
      xray.scatterer( #37
        label='H9B'),
      xray.scatterer( #38
        label='C10',
        site=(-0.444628, 1.167026, 0.565148),
        u=(0.000289, 0.000234, 0.000210,
           0.000018, 0.000036, 0.000040)),
      xray.scatterer( #39
        label='H10B'),
      xray.scatterer( #40
        label='C11',
        site=(-0.371995, 1.272625, 0.676806),
        u=(0.000324, 0.000199, 0.000264,
           0.000054, 0.000089, 0.000009)),
      xray.scatterer( #41
        label='H11'),
      xray.scatterer( #42
        label='C12',
        site=(-0.453277, 1.252159, 0.788677),
        u=(0.000396, 0.000422, 0.000263,
           0.000043, 0.000125, -0.000061)),
      xray.scatterer( #43
        label='H12B'),
      xray.scatterer( #44
        label='H12A')
    )))
  for sc in xs.scatterers():
    sc.flags.set_grad_site(True)
    if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
    if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)

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
        constrained_site_indices=(26, 22)),
      _.tertiary_ch_site(
        pivot=23,
        constrained_site_indices=(24,)),
      _.tertiary_ch_site(
        pivot=25,
        constrained_site_indices=(27,)),
      _.tertiary_ch_site(
        pivot=28,
        constrained_site_indices=(29,)),
      _.tertiary_ch_site(
        pivot=30,
        constrained_site_indices=(31,)),
      _.secondary_ch2_sites(
        pivot=33,
        constrained_site_indices=(35, 34)),
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

  if 0:
    from crys3d.qttbx.xray_structure_viewer import display
    display(xray_structure=xs)

  warned_once = False
  for sc, params in itertools.izip(reparametrisation.structure.scatterers(),
                                   reparametrisation.asu_scatterer_parameters):
    if sc.scattering_type != 'H':
      assert isinstance(params.site, core.independent_site_parameter)
      assert params.site.scatterers[0].label == sc.label
      assert isinstance(params.u, core.independent_cartesian_adp)
      assert params.u.scatterers[0].label == sc.label
      assert isinstance(params.occupancy, core.independent_occupancy_parameter)
      assert params.occupancy.scatterers[0].label == sc.label
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
      assert ([ sc1.label for sc1 in params.site.argument(0).scatterers ]
              == [expected_pivot]), sc.label

  target = (
    59,60,61 , # O1.site
    80,81,82,83,84,85 , # O1.u
    0,1,2 , # O2.site
    86,87,88,89,90,91 , # O2.u
    304,305,306 , # H2.site
    92 , # H2.u
    7,8,9 , # O3.site
    93,94,95,96,97,98 , # O3.u
    307,308,309 , # H3.site
    99 , # H3.u
    14,15,16 , # O4.site
    100,101,102,103,104,105 , # O4.u
    310,311,312 , # H4.site
    106 , # H4.u
    21,22,23 , # O5.site
    107,108,109,110,111,112 , # O5.u
    313,314,315 , # H5.site
    113 , # H5.u
    76,77,78 , # O6.site
    114,115,116,117,118,119 , # O6.u
    28,29,30 , # O7.site
    120,121,122,123,124,125 , # O7.u
    316,317,318 , # H7.site
    126 , # H7.u
    35,36,37 , # O8.site
    127,128,129,130,131,132 , # O8.u
    319,320,321 , # H8.site
    133 , # H8.u
    42,43,44 , # O9.site
    134,135,136,137,138,139 , # O9.u
    322,323,324 , # H9.site
    140 , # H9.u
    49,50,51 , # O10.site
    141,142,143,144,145,146 , # O10.u
    325,326,327 , # H10.site
    147 , # H10.u
    66,67,68 , # O11.site
    148,149,150,151,152,153 , # O11.u
    56,57,58 , # C1.site
    154,155,156,157,158,159 , # C1.u
    328,329,330 , # H1.site
    160 , # H1.u
    3,4,5 , # C2.site
    161,162,163,164,165,166 , # C2.u
    334,335,336 , # H2B.site
    167 , # H2B.u
    10,11,12 , # C3.site
    168,169,170,171,172,173 , # C3.u
    337,338,339 , # H3B.site
    174 , # H3B.u
    17,18,19 , # C4.site
    175,176,177,178,179,180 , # C4.u
    331,332,333 , # H2A.site
    181 , # H2A.u
    340,341,342 , # H4B.site
    182 , # H4B.u
    24,25,26 , # C5.site
    183,184,185,186,187,188 , # C5.u
    343,344,345 , # H5B.site
    189 , # H5B.u
    63,64,65 , # C6.site
    190,191,192,193,194,195 , # C6.u
    346,347,348 , # H6.site
    196 , # H6.u
    69,70,71 , # C7.site
    197,198,199,200,201,202 , # C7.u
    52,53,54 , # C8.site
    203,204,205,206,207,208 , # C8.u
    352,353,354 , # H8A.site
    209 , # H8A.u
    349,350,351 , # H8B.site
    210 , # H8B.u
    45,46,47 , # C9.site
    211,212,213,214,215,216 , # C9.u
    355,356,357 , # H9B.site
    217 , # H9B.u
    38,39,40 , # C10.site
    218,219,220,221,222,223 , # C10.u
    358,359,360 , # H10B.site
    224 , # H10B.u
    73,74,75 , # C11.site
    225,226,227,228,229,230 , # C11.u
    361,362,363 , # H11.site
    231 , # H11.u
    31,32,33 , # C12.site
    232,233,234,235,236,237 , # C12.u
    364,365,366 , # H12B.site
    238 , # H12B.u
    367,368,369 , # H12A.site
    239 , # H12A.u
    )
  assert tuple(reparametrisation.mapping_to_grad_fc) == target,\
         str(reparametrisation)



def run():
  import libtbx.utils
  libtbx.utils.show_times_at_exit()
  exercise()

if __name__ == '__main__':
  run()
