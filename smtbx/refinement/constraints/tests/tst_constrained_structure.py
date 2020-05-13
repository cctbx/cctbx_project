from __future__ import absolute_import, division, print_function
from operator import itemgetter
from cctbx import crystal, xray
from cctbx.array_family import flex
from smtbx.refinement import constraints
import smtbx.refinement.constraints.all as _
import smtbx.refinement.constraints as core
import smtbx.utils
from smtbx.refinement import least_squares
from smtbx import development
from scitbx.lstbx import normal_eqns_solving
from scitbx import matrix
from six.moves import zip
from smtbx.refinement.constraints import fpfdp
from libtbx.test_utils import approx_equal


class test_case(object):

  expected_reparametrisation_for_hydrogen_named = None
  expected_mapping_to_grad_fc = None
  refinement_config = "default"

  def __init__(self, normal_eqns_solving_method):
    self.normal_eqns_solving_method = normal_eqns_solving_method

  def check_reparametrisation_construction(self):
    warned_once = False
    for sc, params in zip(
      self.reparametrisation.structure.scatterers(),
      self.reparametrisation.asu_scatterer_parameters
      ):
      if sc.scattering_type != 'H':
        assert (
          (isinstance(params.site, core.independent_site_parameter)
           and
           isinstance(params.u, core.independent_u_star_parameter)
           )
          or
          (isinstance(params.site, core.special_position_site_parameter)
           and
           isinstance(params.u, core.special_position_u_star_parameter))
          )
        assert params.site.scatterers[0].label == sc.label
        assert params.u.scatterers[0].label == sc.label
        assert isinstance(params.occupancy,
                          core.independent_occupancy_parameter)
        assert params.occupancy.scatterers[0].label == sc.label
      else:
        try:
          (expected_type, expected_pivot) = \
           self.expected_reparametrisation_for_hydrogen_named[sc.label]
          assert isinstance(params.site, expected_type), \
                 (sc.label, params.site, expected_type)
          assert ([ sc1.label for sc1 in params.site.argument(0).scatterers ]
                  == [expected_pivot]), sc.label
        except KeyError:
          if not warned_once:
            print("Warning: incomplete test coverage for H constraint types")
            warned_once = True
            continue
    self.check_reparametrisation_construction_more()

  def check_reparametrisation_construction_more(self):
    """ To be overriden by heirs that needs to perform extra tests """

  def check_mapping_to_grad_fc(self):
    if self.expected_mapping_to_grad_fc is not None:
      assert (tuple(self.reparametrisation.mapping_to_grad_fc)
              == self.expected_mapping_to_grad_fc)
    else:
      print("No mapping to grad Fc test")


  def check_refinement_stability(self):
    if not self.shall_refine_thermal_displacements:
      for sc in self.xray_structure.scatterers():
        sc.flags.set_grad_site(True)
        if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(False)
        if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(False)

    xs = self.xray_structure
    xs0 = self.reference_xray_structure = xs.deep_copy_scatterers()
    mi = xs0.build_miller_set(anomalous_flag=False, d_min=0.5)
    fo_sq = mi.structure_factors_from_scatterers(
      xs0, algorithm="direct").f_calc().norm()
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))

    xs.shake_sites_in_place(rms_difference=0.1)
    if self.shall_refine_thermal_displacements:
      # a spread of 10 for u_iso's would be enormous for our low temperature
      # test structures if those u_iso's were not constrained
      xs.shake_adp(spread=10, # absolute
                   aniso_spread=0.2) # relative

    self.reparametrisation = constraints.reparametrisation(
      xs, self.constraints, self.connectivity_table,
      temperature=self.t_celsius)
    obs = fo_sq.as_xray_observations()
    ls = least_squares.crystallographic_ls(
      obs,
      self.reparametrisation,
      weighting_scheme=least_squares.mainstream_shelx_weighting())
    self.cycles = self.normal_eqns_solving_method(ls)
    print ("%i %s iterations to recover from shaking"
           % (self.cycles.n_iterations,
              self.cycles))
    if 0:
      from crys3d.qttbx.xray_structure_viewer import display
      display(xray_structure=xs)

    diff = xray.meaningful_site_cart_differences(xs0, xs)
    assert diff.max_absolute() < self.site_refinement_tolerance,\
           self.__class__.__name__

    if self.shall_refine_thermal_displacements:
      delta_u = []
      for sc, sc0 in zip(xs.scatterers(), xs0.scatterers()):
        if not sc.flags.use_u_aniso() or not sc0.flags.use_u_aniso(): continue
        delta_u.extend(matrix.col(sc.u_star) - matrix.col(sc0.u_star))
      delta_u = flex.double(delta_u)

      assert flex.max_absolute(delta_u) < self.u_star_refinement_tolerance,\
             self.__class__.__name__


  def display_structure(self):
    from crys3d.qttbx.xray_structure_viewer import display
    display(xray_structure=self.xray_structure)

  def run(self):
    print("[ %s ]" % self.__class__.__name__)
    self.connectivity_table = smtbx.utils.connectivity_table(
      self.xray_structure)
    if self.refinement_config == "default":
      for sc in self.xray_structure.scatterers():
        sc.flags.set_grad_site(True)
        if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
        if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)
    elif self.refinement_config == "fpfdp":
      self.xray_structure.set_inelastic_form_factors(1.54184, 'sasaki')
      for sc in self.xray_structure.scatterers():
        sc.flags.set_grad_site(False)
        sc.flags.set_grad_u_iso(False)
        sc.flags.set_grad_u_aniso(False)
        sc.flags.set_grad_occupancy(False)
        if sc.element_symbol()=='O':
          sc.flags.set_use_fp_fdp(True)
          sc.flags.set_grad_fp(True)
          sc.flags.set_grad_fdp(True)
    elif self.refinement_config == "adp":
      for sc in self.xray_structure.scatterers():
        sc.flags.set_grad_site(False)
        if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
        if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)


    self.reparametrisation = constraints.reparametrisation(
      self.xray_structure,
      self.constraints,
      self.connectivity_table,
      temperature=self.t_celsius,
    )

    self.check_reparametrisation_construction()
    self.check_mapping_to_grad_fc()

    # above settings are customised in the following tests
    # n.b. the fp,fdp tests override the default check_refinement_stability
    self.check_refinement_stability()

class scalar_scaled_adp_test_case(test_case):

  refinement_config = "adp"

  def __init__(self, m):
    test_case.__init__(self, m)
    self.xray_structure = development.sucrose()
    self.t_celsius = 20
    self.shall_refine_thermal_displacements = True
    self.constraints = [
        _.scalar_scaled_u(range(len(self.xray_structure.scatterers())))
        ]
    self.expected_reparametrisation_for_hydrogen_named = {}
    xs = self.xray_structure

    # Make one O atom isotropic to exercise scalar_scaled_u_iso on a non-H atom
    sc1 = xs.scatterers()[1]
    sc1.u_iso = sc1.u_iso_or_equiv(xs.unit_cell())
    sc1.set_use_u_iso_only()
    self.expected_mapping_to_grad_fc = tuple(range(271, 271+155)) # There are
    # 155 adp params: 12 aniso C, 22 isot H, 10 aniso O, 1 isot O

  def check_reparametrisation_construction(self):
    for sc, params in zip(
      self.reparametrisation.structure.scatterers(),
      self.reparametrisation.asu_scatterer_parameters
      ):
      if sc.scattering_type != 'H':
        assert (
            isinstance(params.site, core.independent_site_parameter)
            or
            isinstance(params.site, core.special_position_u_star_parameter)
            )
        if sc.flags.use_u_iso():
          assert isinstance(params.u, core.scalar_scaled_u_iso)
        else:
          assert sc.flags.use_u_aniso()
          assert isinstance(params.u, core.scalar_scaled_u_star)
        assert params.site.scatterers[0].label == sc.label
        assert params.u.scatterers[0].label == sc.label
        assert isinstance(params.occupancy,
                          core.independent_occupancy_parameter)
        assert params.occupancy.scatterers[0].label == sc.label
      else:
        assert isinstance(params.u, core.scalar_scaled_u_iso)
    self.check_reparametrisation_construction_more()

  def check_refinement_stability(self):
    xs = self.xray_structure
    xs0 = self.reference_xray_structure = xs.deep_copy_scatterers()

    # First we construct the Fo array
    mi = xs0.build_miller_set(anomalous_flag=True, d_min=0.5)
    fo_sq = mi.structure_factors_from_scatterers(
        xs0, algorithm="direct").f_calc().norm()
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))

    # Then shake adps and store the shaken adp components
    self.shake_selection = flex.bool([
        sc.element_symbol() in ['C', 'O']
        for sc in xs.scatterers()
        ])
    xs.shake_adp(selection=self.shake_selection)
    adp_ref = []
    for sc in xs.scatterers():
      if sc.flags.use_u_aniso():
        adp_ref.extend(sc.u_star)
      else:
        adp_ref.append(sc.u_iso)

    # Then do the refinement
    self.reparametrisation = constraints.reparametrisation(
      xs, self.constraints, self.connectivity_table,
      temperature=self.t_celsius)
    self.obs = fo_sq.as_xray_observations()
    ls = least_squares.crystallographic_ls(
      self.obs,
      self.reparametrisation,
      weighting_scheme=least_squares.mainstream_shelx_weighting())
    self.cycles = self.normal_eqns_solving_method(ls)
    print ("%i %s iterations to recover from shaking"
           % (self.cycles.n_iterations,
              self.cycles))

    # Then verify the final ADPs have all changed by the same amount
    adp_final = []
    for sc in xs.scatterers():
      if sc.flags.use_u_aniso():
        adp_final.extend(sc.u_star)
      else:
        adp_final.append(sc.u_iso)
    adp_ratios = [x/y for x,y in zip(adp_ref, adp_final) if x>1e-6 and y>1e-6]
    assert all(
        [approx_equal(x, adp_ratios[0], eps=1e-9) for x in adp_ratios[1:]]
    )

    #Make sure the ADPs have changed by some non-zero amount
    assert adp_ratios[0] > 1.01 or adp_ratios[0] < 0.99


class fpfdp_test_case(test_case):

  refinement_config = "fpfdp"

  def __init__(self, m):
    test_case.__init__(self, m)
    self.xray_structure = development.sucrose()

    self.t_celsius = 20
    self.shall_refine_thermal_displacements = False
    self.constraints = []

    self.expected_mapping_to_grad_fc = (
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
        20, 21 # that was easy
        )
    self.expected_reparametrisation_for_hydrogen_named = {}

    self.fp_refinement_tolerance = 1e-6
    self.fdp_refinement_tolerance = 1e-6


  def check_refinement_stability(self):
    xs = self.xray_structure
    xs0 = self.reference_xray_structure = xs.deep_copy_scatterers()
    mi = xs0.build_miller_set(anomalous_flag=True, d_min=0.5)
    fo_sq = mi.structure_factors_from_scatterers(
      xs0, algorithm="direct").f_calc().norm()
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))

    self.shake_selection = flex.bool([
        True if sc.element_symbol()=='O' else False
        for sc
        in self.xray_structure.scatterers()
        ])
    xs.shake_fps(selection=self.shake_selection)
    xs.shake_fdps(selection=self.shake_selection)

    self.reparametrisation = constraints.reparametrisation(
      xs, self.constraints, self.connectivity_table,
      temperature=self.t_celsius)
    self.obs = fo_sq.as_xray_observations()
    ls = least_squares.crystallographic_ls(
      self.obs,
      self.reparametrisation,
      weighting_scheme=least_squares.mainstream_shelx_weighting())
    self.cycles = self.normal_eqns_solving_method(ls)
    print ("%i %s iterations to recover from shaking"
           % (self.cycles.n_iterations,
              self.cycles))

    delta_fp = flex.double([
      sc.fp - sc0.fp
      for sc, sc0
      in zip(xs.scatterers(), xs0.scatterers())
      ])

    delta_fdp = flex.double([
      sc.fdp - sc0.fdp
      for sc, sc0
      in zip(xs.scatterers(), xs0.scatterers())
      ])

    assert flex.max_absolute(delta_fp) < self.fp_refinement_tolerance,\
           self.__class__.__name__

    assert flex.max_absolute(delta_fdp) < self.fdp_refinement_tolerance,\
           self.__class__.__name__


class constrained_fpfdp_test_case(fpfdp_test_case):

  def __init__(self, m):
    fpfdp_test_case.__init__(self, m)
    self.constraints = [
        fpfdp.shared_fp([0,1,3,5,7,9,10,12,14,16,18]),
        fpfdp.shared_fdp([0,1,3,5,7,9,10,12,14,16,18]),
        ]
    self.expected_mapping_to_grad_fc = (
        0, 1, 410, 420, 411, 421, 412, 422, 413, 423, 414, 424, 415, 425, 416,
        426, 417, 427, 418, 428, 419, 429
        )

  def check_reparametrisation_construction_more(self):
    assert self.reparametrisation.n_independents==2

  def check_refinement_stability(self):
    """Shake the fps and fdps, then refine, then make sure they are all equal
    within a very small tolerance.
    """
    fpfdp_test_case.check_refinement_stability(self)
    oxygen_fps = [
        sc.fp for sc in self.xray_structure.scatterers()
        if sc.element_symbol()=='O'
        ]
    oxygen_fdps = [
        sc.fdp for sc in self.xray_structure.scatterers()
        if sc.element_symbol()=='O'
        ]
    assert len([
      x for x in oxygen_fps
      if approx_equal(x, oxygen_fps[0], eps=1e-9, multiplier=None)
      ]) == len(oxygen_fps)
    assert len([
      x for x in oxygen_fdps
      if approx_equal(x, oxygen_fdps[0], eps=1e-9, multiplier=None)
      ]) == len(oxygen_fdps)


class sucrose_test_case(test_case):
  """
  sucrose from Olex 2 samples

  Notes:
  - atom H2A has been moved down the list to test non-contiguous indices
    in argument 'constrained_site_indices'

  - the sites of H12A and H12B have been swapped (this is purely conventional:
    it turns out that the smtbx code does not follow that of ShelXL which
    was used to produce this structure in the first place)
  """

  def __init__(self, m):
    test_case.__init__(self, m)
    self.xray_structure = development.sucrose()

    self.t_celsius = 20
    self.shall_refine_thermal_displacements = False

    self.constraints = [
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=1,
        constrained_site_indices=(2,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=3,
        constrained_site_indices=(4,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=5,
        constrained_site_indices=(6,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=7,
        constrained_site_indices=(8,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=10,
        constrained_site_indices=(11,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=12,
        constrained_site_indices=(13,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=14,
        constrained_site_indices=(15,)),
      _.terminal_tetrahedral_xh_site(
        rotating=True,
        pivot=16,
        constrained_site_indices=(17,)),
      _.tertiary_xh_site(
        pivot=19,
        constrained_site_indices=(20,)),
      _.secondary_xh2_sites(
        pivot=21,
        flapping = True,
        constrained_site_indices=(26, 22)),
      _.tertiary_xh_site(
        pivot=23,
        constrained_site_indices=(24,)),
      _.tertiary_xh_site(
        pivot=25,
        constrained_site_indices=(27,)),
      _.tertiary_xh_site(
        pivot=28,
        constrained_site_indices=(29,)),
      _.tertiary_xh_site(
        pivot=30,
        constrained_site_indices=(31,)),
      _.secondary_xh2_sites(
        pivot=33,
        flapping = True,
        constrained_site_indices=(35, 34)),
      _.tertiary_xh_site(
        pivot=36,
        constrained_site_indices=(37,)),
      _.tertiary_xh_site(
        pivot=38,
        constrained_site_indices=(39,)),
      _.tertiary_xh_site(
        pivot=40,
        constrained_site_indices=(41,)),
      _.secondary_xh2_sites(
        pivot=42,
        flapping = True,
        constrained_site_indices=(43, 44)),
      ]

    self.expected_reparametrisation_for_hydrogen_named = {
      "H2": (core.terminal_tetrahedral_xh_site, 'O2'),
      "H3": (core.terminal_tetrahedral_xh_site, 'O3'),
      "H4": (core.terminal_tetrahedral_xh_site, 'O4'),
      "H5": (core.terminal_tetrahedral_xh_site, 'O5'),
      "H7": (core.terminal_tetrahedral_xh_site, 'O7'),
      "H8": (core.terminal_tetrahedral_xh_site, 'O8'),
      "H9": (core.terminal_tetrahedral_xh_site, 'O9'),
      "H10": (core.terminal_tetrahedral_xh_site, 'O10'),
      "H1": (core.tertiary_xh_site, 'C1'),
      "H2A": (core.secondary_xh2_sites, 'C2'),
      "H2B": (core.secondary_xh2_sites, 'C2'),
      "H3B": (core.tertiary_xh_site, 'C3'),
      "H4B": (core.tertiary_xh_site, 'C4'),
      "H5B": (core.tertiary_xh_site, 'C5'),
      "H6": (core.tertiary_xh_site, 'C6'),
      "H8A": (core.secondary_xh2_sites, 'C8'),
      "H8B": (core.secondary_xh2_sites, 'C8'),
      "H9B": (core.tertiary_xh_site, 'C9'),
      "H10B": (core.tertiary_xh_site, 'C10'),
      "H11": (core.tertiary_xh_site, 'C11'),
      "H12A": (core.secondary_xh2_sites, 'C12'),
      "H12B": (core.secondary_xh2_sites, 'C12'),
    }

    self.expected_mapping_to_grad_fc = (
            59, 60, 61, # O1.site
            80, 81, 82, 83, 84, 85, # O1.u
            0, 1, 2, # O2.site
            86, 87, 88, 89, 90, 91,# O2.u
            394, 395, 396, # H2.site
            92, # H2.u
            7, 8, 9, # O3.site
            93, 94, 95, 96, 97, 98, # O3.u
            397, 398, 399, # H3.site
            99, # H3.u
            14, 15, 16, # O4.site
            100, 101, 102, 103, 104, 105, # O4.u
            400, 401, 402, # H4.site
            106, # H4.u
            21, 22, 23, # O5.site
            107, 108, 109, 110, 111, 112, # O5.u
            403, 404, 405, # H5.site
            113, # H5.u
            76, 77, 78, # O6.site
            114, 115, 116, 117, 118, 119, # O6.u
            28, 29, 30, # O7.site
            120, 121, 122, 123, 124, 125, # O7.u
            406, 407, 408, # H7.site
            126, # H7.u
            35, 36, 37, # O8.site
            127, 128, 129, 130, 131, 132, # O8.u
            409, 410, 411, # H8.site
            133, # H8.u
            42, 43, 44, # O9.site
            134, 135, 136, 137, 138, 139, # O9.u
            412, 413, 414, # H9.site
            140, # H9.u
            49, 50, 51, # O10.site
            141, 142, 143, 144, 145, 146, # O10.u
            415, 416, 417, # H10.site
            147, # H10.u
            66, 67, 68, # O11.site
            148, 149, 150, 151, 152, 153, # O11.u
            56, 57, 58, # C1.site
            154, 155, 156, 157, 158, 159, # C1.u
            418, 419, 420, # H1.site
            160, # H1.u
            3, 4, 5, # C2.site
            161, 162, 163, 164, 165, 166, # C2.u
            424, 425, 426, # H2B.site
            167, # H2B.u
            10, 11, 12, # C3.site
            168, 169, 170, 171, 172, 173, # C3.u
            427, 428, 429, # H3B.site
            174, # H3B.u
            17, 18, 19, # C4.site
            175, 176, 177, 178, 179, 180, # C4.u
            421, 422, 423, # H2A.site
            181, # H2A.u
            430, 431, 432, # H4B.site
            182, # H4B.u
            24, 25, 26, # C5.site
            183, 184, 185, 186, 187, 188, # C5.u
            433, 434, 435, # H5B.site
            189, # H5B.u
            63, 64, 65, # C6.site
            190, 191, 192, 193, 194, 195, # C6.u
            436, 437, 438, # H6.site
            196, # H6.u
            69, 70, 71, # C7.site
            197, 198, 199, 200, 201, 202, # C7.u
            52, 53, 54, # C8.site
            203, 204, 205, 206, 207, 208, # C8.u
            442, 443, 444, # H8A.site
            209, # H8A.u
            439, 440, 441, # H8B.site
            210, # H8B.u
            45, 46, 47, # C9.site
            211, 212, 213, 214, 215, 216, # C9.u
            445, 446, 447, # H9B.site
            217, # H9B.u
            38, 39, 40, # C10.site
            218, 219, 220, 221, 222, 223, # C10.u
            448, 449, 450, # H10B.site
            224, # H10B.u
            73, 74, 75, # C11.site
            225, 226, 227, 228, 229, 230, # C11.u
            451, 452, 453, # H11.site
            231, # H11.u
            31, 32, 33, # C12.site
            232, 233, 234, 235, 236, 237, # C12.u
            454, 455, 456, # H12B.site
            238, # H12B.u
            457, 458, 459, # H12A.site
            239 # H12A.u
            )

    self.site_refinement_tolerance = 1e-4


class saturated_test_case(test_case):
  """ Durham database: 03srv020

  H1N and H2N sites and u's have been swapped.
  """

  def __init__(self, m):
    test_case.__init__(self, m)
    self.xray_structure = xray.structure(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(3.753, 14.54, 15.868, 90, 92.58, 90),
        space_group_symbol='hall: -P 2ybc (x-z,y,z)'),
      scatterers=flex.xray_scatterer((
        xray.scatterer( #0
                        label='O1',
                        site=(0.299733, 0.262703, 0.397094),
                        u=(0.003622, 0.000123, 0.000108,
                           0.000154, -0.000304, -0.000017)),
        xray.scatterer( #1
                        label='O2',
                        site=(0.606432, 0.145132, 0.437285),
                        u=(0.004117, 0.000149, 0.000118,
                           0.000405, -0.000117, -0.000042)),
        xray.scatterer( #2
                        label='N1',
                        site=(0.481175, 0.221358, 0.451529),
                        u=(0.001750, 0.000091, 0.000090,
                           0.000016, -0.000027, -0.000000)),
        xray.scatterer( #3
                        label='N2',
                        site=(0.669716, 0.393801, 0.763893),
                        u=(0.002874, 0.000122, 0.000071,
                           0.000168, -0.000071, -0.000000)),
        xray.scatterer( #4
                        label='H1N',
                        site=(0.777763, 0.365311, 0.806784),
                        u=0.042273),
        xray.scatterer( #5
                        label='H2N',
                        site=(0.589373, 0.450143, 0.770096),
                        u=0.042273),
        xray.scatterer( #6
                        label='C1',
                        site=(0.542801, 0.263225, 0.532794),
                        u=(0.001423, 0.000084, 0.000076,
                           -0.000009, 0.000000, -0.000003)),
        xray.scatterer( #7
                        label='C2',
                        site=(0.718203, 0.216166, 0.600719),
                        u=(0.001241, 0.000081, 0.000088,
                           -0.000011, 0.000023, 0.000009)),
        xray.scatterer( #8
                        label='C3',
                        site=(0.754550, 0.260840, 0.677732),
                        u=(0.001553, 0.000097, 0.000079,
                           0.000044, -0.000014, 0.000019)),
        xray.scatterer( #9
                        label='H3',
                        site=(0.868083, 0.229804, 0.724284),
                        u=0.031215),
        xray.scatterer( #10
                        label='C4',
                        site=(0.627437, 0.351296, 0.688971),
                        u=(0.001601, 0.000093, 0.000075,
                           0.000015, 0.000027, 0.000004)),
        xray.scatterer( #11
                        label='C5',
                        site=(0.454744, 0.396843, 0.619275),
                        u=(0.001302, 0.000082, 0.000083,
                           -0.000002, -0.000011, 0.000007)),
        xray.scatterer( #12
                        label='C6',
                        site=(0.414545, 0.352153, 0.542558),
                        u=(0.001310, 0.000086, 0.000080,
                           0.000014, -0.000038, 0.000012)),
        xray.scatterer( #13
                        label='H6',
                        site=(0.298067, 0.382503, 0.496003),
                        u=0.028457),
        xray.scatterer( #14
                        label='C7',
                        site=(0.870740, 0.125780, 0.594964),
                        u=(0.001485, 0.000100, 0.000087,
                           0.000015, 0.000006, 0.000009)),
        xray.scatterer( #15
                        label='C8',
                        site=(1.017756, 0.053520, 0.594148),
                        u=(0.002161, 0.000095, 0.000112,
                           0.000070, 0.000030, 0.000015)),
        xray.scatterer( #16
                        label='H8',
                        site=(1.135411, -0.004309, 0.593494),
                        u=0.039284),
        xray.scatterer( #17
                        label='C9',
                        site=(0.317957, 0.487839, 0.631013),
                        u=(0.001562, 0.000101, 0.000080,
                           0.000008, -0.000022, 0.000009)),
        xray.scatterer( #18
                        label='C10',
                        site=(0.204108, 0.561746, 0.646136),
                        u=(0.001931, 0.000100, 0.000123,
                           0.000052, -0.000025, -0.000002)),
        xray.scatterer( #19
                        label='H10',
                        site=(0.112835, 0.620998, 0.658260),
                        u=0.039775)
      )))


    self.t_celsius = -153
    self.shall_refine_thermal_displacements = True

    k=1.5 # that is the multiplier used to refine the structure with ShelXL
    self.constraints = [
      _.terminal_planar_xh2_sites(
        pivot=3,
        constrained_site_indices=(4, 5)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=3,
        u_iso_scatterer_idx=4,
        multiplier=k),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=3,
        u_iso_scatterer_idx=5,
        multiplier=k),

      _.terminal_linear_ch_site(
        pivot=15,
        constrained_site_indices=(16,)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=15,
        u_iso_scatterer_idx=16,
        multiplier=k),

      _.terminal_linear_ch_site(
        pivot=18,
        constrained_site_indices=(19,)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=18,
        u_iso_scatterer_idx=19,
        multiplier=k),

      _.secondary_planar_xh_site(
        pivot=12,
        constrained_site_indices=(13,)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=12,
        u_iso_scatterer_idx=13,
        multiplier=k),

      _.secondary_planar_xh_site(
        pivot=8,
        constrained_site_indices=(9,)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=8,
        u_iso_scatterer_idx=9,
        multiplier=k),
      ]

    self.expected_reparametrisation_for_hydrogen_named = {
      'H1N': (core.terminal_planar_xh2_sites, 'N2'),
      'H2N': (core.terminal_planar_xh2_sites, 'N2'),
      'H10': (core.terminal_linear_ch_site, 'C10'),
      'H6' : (core.secondary_planar_xh_site, 'C6'),
      'H3' : (core.secondary_planar_xh_site, 'C3'),
      'H8': (core.terminal_linear_ch_site, 'C8'),
    }

    self.site_refinement_tolerance = 1e-2
    self.u_star_refinement_tolerance = 1e-5


class symmetry_equivalent_test_case(test_case):
  """ 09srv172 from Durham database """

  def __init__(self, m):
    test_case.__init__(self, m)
    self.xray_structure = xray.structure(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(17.0216, 8.4362, 10.2248, 90, 102.79, 90),
        space_group_symbol='hall: -C 2yc'),
      scatterers=flex.xray_scatterer((
        xray.scatterer( #0
                        label='S1',
                        site=(0.525736, 0.737492, 0.619814),
                        u=(0.000084, 0.000243, 0.000191,
                           0.000021, 0.000026, -0.000035)),
        xray.scatterer( #1
                        label='C1',
                        site=(0.500000, 0.868009, 0.750000),
                        u=(0.000061, 0.000181, 0.000164,
                           0.000000, 0.000025, 0.000000)),
        xray.scatterer( #2
                        label='C2',
                        site=(0.533017, 0.552825, 0.710913),
                        u=(0.000161, 0.000241, 0.000348,
                           0.000041, 0.000039, -0.000020)),
        xray.scatterer( #3
                        label='H2A',
                        site=(0.525986, 0.462184, 0.647894),
                        u=0.041420),
        xray.scatterer( #4
                        label='H2B',
                        site=(0.586473, 0.543642, 0.772877),
                        u=0.038360),
        xray.scatterer( #5
                        label='C3',
                        site=(0.425914, 0.971290, 0.682589),
                        u=(0.000058, 0.000199, 0.000164,
                           -0.000003, 0.000017, -0.000007)),
        xray.scatterer( #6
                        label='H3',
                        site=(0.441258, 1.029902, 0.606966),
                        u=0.023950),
        xray.scatterer( #7
                        label='C4',
                        site=(0.349971, 0.874741, 0.622481),
                        u=(0.000064, 0.000236, 0.000219,
                           -0.000015, 0.000011, -0.000014)),
        xray.scatterer( #8
                        label='H4B',
                        site=(0.362228, 0.799281, 0.555566),
                        u=0.026970),
        xray.scatterer( #9
                        label='H4A',
                        site=(0.333906, 0.812566, 0.694426),
                        u=0.025070),
        xray.scatterer( #10
                        label='C5',
                        site=(0.279832, 0.981636, 0.555089),
                        u=(0.000064, 0.000307, 0.000222,
                           -0.000007, 0.000003, -0.000021)),
        xray.scatterer( #11
                        label='H5B',
                        site=(0.294150, 1.037327, 0.478372),
                        u=0.026970),
        xray.scatterer( #12
                        label='H5A',
                        site=(0.231706, 0.915594, 0.519984),
                        u=0.034720),
        xray.scatterer( #13
                        label='C6',
                        site=(0.259978, 1.103478, 0.653453),
                        u=(0.000061, 0.000335, 0.000246,
                           0.000016, 0.000021, -0.000001)),
        xray.scatterer( #14
                        label='H6B',
                        site=(0.216403, 1.174193, 0.606307),
                        u=0.037560),
        xray.scatterer( #15
                        label='H6A',
                        site=(0.240791, 1.048477, 0.726038),
                        u=0.024860),
        xray.scatterer( #16
                        label='C7',
                        site=(0.334517, 1.201490, 0.713404),
                        u=(0.000071, 0.000265, 0.000251,
                           0.000019, 0.000019, -0.000029)),
        xray.scatterer( #17
                        label='H7B',
                        site=(0.349843, 1.265431, 0.641770),
                        u=0.031500),
        xray.scatterer( #18
                        label='H7A',
                        site=(0.321733, 1.275376, 0.780996),
                        u=0.030620),
        xray.scatterer( #19
                        label='C8',
                        site=(0.405687, 1.096392, 0.779637),
                        u=(0.000064, 0.000247, 0.000187,
                           0.000009, 0.000015, -0.000033)),
        xray.scatterer( #20
                        label='H8A',
                        site=(0.392638, 1.042204, 0.858057),
                        u=0.026000),
        xray.scatterer( #21
                        label='H8B',
                        site=(0.453574, 1.163891, 0.812432),
                        u=0.027360)
      )))

    self.t_celsius = -153
    self.shall_refine_thermal_displacements = True

    self.constraints = [
      _.secondary_xh2_sites(
        pivot=2,
        constrained_site_indices=(3,4)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=2,
        u_iso_scatterer_idx=3,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=2,
        u_iso_scatterer_idx=4,
        multiplier=1.5),

      _.tertiary_xh_site(
        pivot=5,
        constrained_site_indices=(6,)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=5,
        u_iso_scatterer_idx=6,
        multiplier=1.5),

      _.secondary_xh2_sites(
        pivot=7,
        constrained_site_indices=(8, 9)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=7,
        u_iso_scatterer_idx=8,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=7,
        u_iso_scatterer_idx=9,
        multiplier=1.5),

      _.secondary_xh2_sites(
        pivot=10,
        constrained_site_indices=(11, 12)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=10,
        u_iso_scatterer_idx=11,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=10,
        u_iso_scatterer_idx=12,
        multiplier=1.5),

      _.secondary_xh2_sites(
        pivot=13,
        constrained_site_indices=(14, 15)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=13,
        u_iso_scatterer_idx=14,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=13,
        u_iso_scatterer_idx=15,
        multiplier=1.5),

      _.secondary_xh2_sites(
        pivot=16,
          constrained_site_indices=(17, 18)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=16,
        u_iso_scatterer_idx=17,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=16,
        u_iso_scatterer_idx=18,
        multiplier=1.5),

      _.secondary_xh2_sites(
        pivot=19,
        constrained_site_indices=(20, 21)),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=19,
        u_iso_scatterer_idx=20,
        multiplier=1.5),

      _.u_iso_proportional_to_pivot_u_eq(
        u_eq_scatterer_idx=19,
        u_iso_scatterer_idx=21,
        multiplier=1.5),
    ]

    self.expected_reparametrisation_for_hydrogen_named = {
      "H2A": (core.secondary_xh2_sites, 'C2'),
      "H2B": (core.secondary_xh2_sites, 'C2'),
      "H3" : (core.tertiary_xh_site   , 'C3'),
      "H4A": (core.secondary_xh2_sites, 'C4'),
      "H4B": (core.secondary_xh2_sites, 'C4'),
      "H5A": (core.secondary_xh2_sites, 'C5'),
      "H5B": (core.secondary_xh2_sites, 'C5'),
      "H6A": (core.secondary_xh2_sites, 'C6'),
      "H6B": (core.secondary_xh2_sites, 'C6'),
      "H7A": (core.secondary_xh2_sites, 'C7'),
      "H7B": (core.secondary_xh2_sites, 'C7'),
      "H8A": (core.secondary_xh2_sites, 'C8'),
      "H8B": (core.secondary_xh2_sites, 'C8'),
    }

    self.site_refinement_tolerance = 0.01
    self.u_star_refinement_tolerance = 5e-7

  def check_reparametrisation_construction_more(self):
    for params in self.reparametrisation.asu_scatterer_parameters:
      if params.site.scatterers[0].label == 'H2A':
        h2a = params.site
        (pivot, pivot_neighbour_0, pivot_neighbour_1,
         bond_length, h_c_h_angle) = h2a.arguments()
        expected = sorted([ (core.independent_site_parameter, 'S1'),
                            (core.symmetry_equivalent_site_parameter, 'C2',
                                '-x+1,y,-z+3/2') ], key=itemgetter(1))
        actual = []
        for n in (pivot_neighbour_0, pivot_neighbour_1):
          if type(n) == core.independent_site_parameter:
            actual.append((type(n), n.scatterers[0].label))
          elif type(n) == core.symmetry_equivalent_site_parameter:
            actual.append((type(n),
                           n.original.scatterers[0].label,
                           str(n.motion)))
        actual.sort(key=itemgetter(1))
        assert actual == expected

def run():
  import libtbx.utils
  libtbx.utils.show_times_at_exit()
  import sys
  from libtbx.option_parser import option_parser
  command_line = (option_parser()
    .option(None, "--normal_eqns_solving_method",
            default='naive')
    .option(None, "--fix_random_seeds",
            action='store_true',
            default='naive')
  ).process(args=sys.argv[1:])
  opts = command_line.options
  if opts.fix_random_seeds:
    import random
    random.seed(1)
    flex.set_random_seed(1)
  gradient_threshold=1e-8
  step_threshold=1e-8
  if opts.normal_eqns_solving_method == 'naive':
    m = lambda eqns: normal_eqns_solving.naive_iterations(
      eqns,
      gradient_threshold=gradient_threshold,
      step_threshold=step_threshold)
  elif opts.normal_eqns_solving_method == 'levenberg-marquardt':
    m = lambda eqns: normal_eqns_solving.levenberg_marquardt_iterations(
      eqns,
      gradient_threshold=gradient_threshold,
      step_threshold=gradient_threshold,
      tau=1e-7)
  else:
    raise RuntimeError("Unknown method %s" % opts.normal_eqns_solving_method)
  for t in [
    saturated_test_case(m),
    sucrose_test_case(m),
    symmetry_equivalent_test_case(m),
    fpfdp_test_case(m),
    constrained_fpfdp_test_case(m),
    scalar_scaled_adp_test_case(m),
    ]:
    t.run()

if __name__ == '__main__':
  run()
