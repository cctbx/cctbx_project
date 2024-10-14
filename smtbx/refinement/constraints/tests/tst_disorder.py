from __future__ import absolute_import, division, print_function

from six.moves import zip
from six.moves import cStringIO as StringIO

disordered = """
TITL 00asb006 in Pca2(1)
CELL 0.71073 25.392 25.853 24.736 90 90 90
ZERR 4 0.008 0.008 0.008 0 0 0
LATT -1
SYMM -X,-Y,0.5+Z
SYMM 0.5+X,-Y,+Z
SYMM 0.5-X,+Y,0.5+Z
SFAC C H N O Co
UNIT 540 920 24 184 24

TEMP -153

WGHT    0.100000    5.070000
FVAR       0.07939

C11   1    0.483531    0.459094    0.598199    11.00000    0.03813    0.02398 =
         0.01913   -0.00385    0.00456   -0.01024
C12   1    0.470584    0.441994    0.539618    11.00000    0.05678    0.04644 =
         0.02077   -0.00864   -0.00707    0.01456
PART 1
C13A  1    0.467614    0.382305    0.543321    10.50000    0.08323
AFIX 137
H13A  2    0.500814    0.368863    0.557836    10.50000   -1.50000
H13B  2    0.461581    0.367872    0.507215    10.50000   -1.50000
H13C  2    0.438564    0.372387    0.567257    10.50000   -1.50000
AFIX   0
C14A  1    0.512981    0.456421    0.502662    10.50000    0.13294
AFIX 137
H14A  2    0.523457    0.492309    0.509528    10.50000   -1.50000
H14B  2    0.500683    0.453086    0.465271    10.50000   -1.50000
H14C  2    0.543247    0.433548    0.508394    10.50000   -1.50000
AFIX   0
C15A  1    0.415133    0.458461    0.522044    10.50000    0.03863
AFIX 137
H15A  2    0.388967    0.441779    0.545370    10.50000   -1.50000
H15B  2    0.409209    0.448013    0.484449    10.50000   -1.50000
H15C  2    0.411795    0.496114    0.525084    10.50000   -1.50000
AFIX   0
PART 0
PART 2
C13B  1    0.499998    0.393168    0.524392    10.50000    0.03074
AFIX 137
H13D  2    0.537495    0.397446    0.532638    10.50000   -1.50000
H13E  2    0.495628    0.386533    0.485644    10.50000   -1.50000
H13F  2    0.485905    0.363926    0.545039    10.50000   -1.50000
AFIX   0
C14B  1    0.496158    0.488076    0.502627    10.50000    0.05672
AFIX 137
H14D  2    0.486815    0.521767    0.518047    10.50000   -1.50000
H14E  2    0.482406    0.485655    0.465695    10.50000   -1.50000
H14F  2    0.534566    0.484402    0.501987    10.50000   -1.50000
AFIX   0
C15B  1    0.416045    0.429952    0.540272    10.50000    0.07696
AFIX 137
H15D  2    0.408861    0.405168    0.569328    10.50000   -1.50000
H15E  2    0.405938    0.414768    0.505484    10.50000   -1.50000
H15F  2    0.395683    0.461606    0.546453    10.50000   -1.50000
AFIX   0
PART 0
HKLF 4
"""

def exercise_with_disorder():
  from iotbx import shelx
  from smtbx import refinement
  import smtbx.refinement.constraints as core

  p = shelx.parse_smtbx_refinement_model(file=StringIO(disordered))
  xs = p.structure
  mi = xs.build_miller_set(d_min=0.6, anomalous_flag=True)
  fcalc = mi.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
  xm = refinement.model(fo_sq=fcalc.norm(),
                        xray_structure=xs,
                        constraints=p.constraints,
                        restraints_manager=p.restraints_manager,
                        weighting_scheme=p.weighting_scheme,
                        temperature_in_celsius=p.temperature_in_celsius,
                        conformer_indices=p.conformer_indices)
  assert xm.temperature_in_celsius == -153
  ls = xm.least_squares()
  expected_reparametrisation_for = {}
  for i in (13, 14, 15):
    for s in 'ABCDEF':
      pivot_letter = 'A' if s < 'D' else 'B'
      expected_reparametrisation_for["H%i%s" % (i,s)] = (
        core.terminal_tetrahedral_xh3_sites, "C%i%s" % (i,pivot_letter))
  for sc, params in zip(
    ls.reparametrisation.structure.scatterers(),
    ls.reparametrisation.asu_scatterer_parameters):
    if sc.scattering_type != 'H':
      continue
    (expected_type, expected_pivot) = expected_reparametrisation_for[sc.label]
    assert isinstance(params.site, expected_type)
    assert ([sc1.label for sc1 in params.site.argument(0).scatterers] ==
            [expected_pivot])

if __name__ == '__main__':
  exercise_with_disorder()
  print('OK')
