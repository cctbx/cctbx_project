from mmtbx import masks
from cctbx.eltbx import van_der_waals_radii
from cctbx import maptbx
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import math
import sys

def structure_init(site,sg,cell):
  symmetry = crystal.symmetry(unit_cell=cell,
                              space_group_symbol=sg)
  structure = xray.structure(crystal_symmetry=symmetry)
  scatterer = xray.scatterer(
                 site = site,
                 u = 0.1,
                 occupancy = 1.0,
                 scattering_type = "C")
  structure.add_scatterer(scatterer)
  xyzf = flex.vec3_double()
  atmrad = flex.double()
  for scatterer in structure.scatterers():
    xyzf.append(list(scatterer.site))
    atmrad.append(van_der_waals_radii.vdw.table[scatterer.element_symbol()])
  assert xyzf.size() == atmrad.size()
  return structure, xyzf, atmrad

def exercise_1():
#------------------------------------------------------------------------TEST-1
 if(1):
  r=(1.67,0.68,0.0,0.5,-0.31,-0.64,-1.63)
  d = 2.0*van_der_waals_radii.vdw.table["C"]
  a = (d/2.) * (4./3.*math.pi*2.)**(1./3.)
  for x in r:
   for y in r:
    for z in r:
     structure, xyzf, atmrad = structure_init((x,y,z),
                                              "P1",
                                              (a, a, a, 90, 90, 90))
     step = 0.1
     crystal_gridding = maptbx.crystal_gridding(
       unit_cell = structure.unit_cell(),
       step = step)
     shrink_truncation_radius = 0.0
     solvent_radius = 0.0
     m = masks.around_atoms(
       structure.unit_cell(),
       structure.space_group().order_z(),
       xyzf,
       atmrad,
       crystal_gridding.n_real(),
       solvent_radius,
       shrink_truncation_radius)
     assert m.solvent_radius == 0
     assert m.shrink_truncation_radius == 0
     assert approx_equal(1./m.accessible_surface_fraction, 2.0, 1e-2)
     assert approx_equal(1./m.contact_surface_fraction, 2.0, 1e-2)
     assert flex.max(m.data) == 1
     assert flex.min(m.data) == 0
     assert m.data.size() == m.data.count(1) + m.data.count(0)
#------------------------------------------------------------------------TEST-2
 if(1):
  r=(-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0)
  asf=[]
  csf=[]
  for x in r:
   for y in r:
     for z in r:
       structure, xyzf, atmrad = structure_init((x,y,z),
                                                "C2",
                                                (5.0, 5.5, 7.0, 90, 80, 90))
       step = 0.25
       crystal_gridding = maptbx.crystal_gridding(
         unit_cell = structure.unit_cell(),
         step = step)
       shrink_truncation_radius = 0.0
       solvent_radius = 0.0
       m = masks.around_atoms(
         structure.unit_cell(),
         structure.space_group().order_z(),
         xyzf,
         atmrad,
         crystal_gridding.n_real(),
         solvent_radius,
         shrink_truncation_radius)
       asf.append(m.accessible_surface_fraction)
       csf.append(m.contact_surface_fraction)
       assert m.accessible_surface_fraction == m.contact_surface_fraction
  assert approx_equal(min(asf),max(asf))
#------------------------------------------------------------------------TEST-3
 if(1):
  r=(1.673,0.685,-0.315,-0.649,-1.637)
  asf=[]
  csf=[]
  for x in r:
   for y in r:
     for z in r:
       structure, xyzf, atmrad = structure_init((x,y,z),
                                                "P1",
                                                (5.0, 5.5, 6.0, 70, 80, 100))
       step = 0.1
       crystal_gridding = maptbx.crystal_gridding(
         unit_cell = structure.unit_cell(),
         step = step)
       shrink_truncation_radius = 0.0
       solvent_radius = 0.0
       m = masks.around_atoms(
          structure.unit_cell(),
          structure.space_group().order_z(),
          xyzf,
          atmrad,
          crystal_gridding.n_real(),
          solvent_radius,
          shrink_truncation_radius)
       asf.append(m.accessible_surface_fraction)
       csf.append(m.contact_surface_fraction)
       assert m.accessible_surface_fraction == m.contact_surface_fraction
  assert approx_equal(min(asf),max(asf), 1e-3)

def exercise_2():
  symmetry = crystal.symmetry(
    unit_cell=(5.67, 10.37, 10.37, 90, 135.49, 90),
    space_group_symbol="C2")
  structure = xray.structure(crystal_symmetry=symmetry)
  atmrad = flex.double()
  xyzf = flex.vec3_double()
  for k in xrange(100):
    scatterer = xray.scatterer(
      site = ((1.+k*abs(math.sin(k)))/1000.0,
              (1.+k*abs(math.cos(k)))/1000.0,
              (1.+ k)/1000.0),
      scattering_type = "C")
    structure.add_scatterer(scatterer)
    atmrad.append(van_der_waals_radii.vdw.table[scatterer.element_symbol()])
    xyzf.append(scatterer.site)
  miller_set = miller.build_set(
    crystal_symmetry=structure,
    d_min=1.0,
    anomalous_flag=False)
  step = 0.5
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell=structure.unit_cell(),
    step=step)
  nxyz = crystal_gridding.n_real()
  shrink_truncation_radius = 1.0
  solvent_radius = 1.0
  m1 = masks.around_atoms(
    structure.unit_cell(),
    structure.space_group().order_z(),
    structure.sites_frac(),
    atmrad,
    nxyz,
    solvent_radius,
    shrink_truncation_radius)
  assert m1.solvent_radius == 1
  assert m1.shrink_truncation_radius == 1
  assert flex.max(m1.data) == 1
  assert flex.min(m1.data) == 0
  assert m1.data.size() == m1.data.count(1) + m1.data.count(0)
  m2 = masks.bulk_solvent(
    xray_structure=structure,
    gridding_n_real=nxyz,
    solvent_radius=solvent_radius,
    shrink_truncation_radius=shrink_truncation_radius)
  assert m2.data.all_eq(m1.data)
  m3 = masks.bulk_solvent(
    xray_structure=structure,
    grid_step=step,
    solvent_radius=solvent_radius,
    shrink_truncation_radius=shrink_truncation_radius)
  assert m3.data.all_eq(m1.data)
  f_mask2 = m2.structure_factors(miller_set=miller_set)
  f_mask3 = m3.structure_factors(miller_set=miller_set)
  assert approx_equal(f_mask2.data(), f_mask3.data())
  assert approx_equal(flex.sum(flex.abs(f_mask3.data())), 1095.17999134)

def exercise_centrics(space_group_info, n_sites=10):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=(("O","N","C")*(n_sites/3+1))[:n_sites],
    volume_per_atom=30,
    min_distance=1)
  for anomalous_flag in [False, True]:
    miller_set = miller.build_set(
      crystal_symmetry=structure,
      d_min=1,
      anomalous_flag=anomalous_flag)
    for shrink_truncation_radius in [0, .5*6**.5]:
      for solvent_radius in [0, .5*5**.5]:
        bulk_solvent_mask = masks.bulk_solvent(
          xray_structure=structure,
          grid_step=0.5,
          solvent_radius=solvent_radius,
          shrink_truncation_radius=shrink_truncation_radius)
        f_mask = bulk_solvent_mask.structure_factors(miller_set=miller_set)
        centrics = f_mask.select_centric()
        if (centrics.indices().size() > 0):
          ideal = centrics.phase_transfer(centrics)
          assert flex.max(flex.abs(ideal.data() - centrics.data())) < 1.e-6

def run_call_back(flags, space_group_info):
  exercise_centrics(space_group_info)

def run():
  exercise_1()
  exercise_2()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
