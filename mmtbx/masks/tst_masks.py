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
    ignore_zero_occupancy_atoms = False,
    solvent_radius=solvent_radius,
    shrink_truncation_radius=shrink_truncation_radius)
  assert m2.data.all_eq(m1.data)
  m3 = masks.bulk_solvent(
    xray_structure=structure,
    grid_step=step,
    ignore_zero_occupancy_atoms = False,
    solvent_radius=solvent_radius,
    shrink_truncation_radius=shrink_truncation_radius)
  assert m3.data.all_eq(m1.data)
  f_mask2 = m2.structure_factors(miller_set=miller_set)
  f_mask3 = m3.structure_factors(miller_set=miller_set)
  assert approx_equal(f_mask2.data(), f_mask3.data())
  assert approx_equal(flex.sum(flex.abs(f_mask3.data())), 1095.17999134)

def exercise_3():
  from mmtbx import masks
  from cctbx import sgtbx
  xs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["C"]*2,
    unit_cell        = (10, 20, 30, 70, 80, 120))
  f_calc = xs.structure_factors(d_min = 2.0).f_calc()
  mp = masks.mask_master_params.extract()
  mm1 = masks.manager(miller_array   = f_calc,
                      xray_structure = xs)
  assert list(xs.scatterers().extract_occupancies()) == [1.0, 1.0]
  fmasks1 = mm1.shell_f_masks()
  assert len(fmasks1) == 1
  assert flex.mean( flex.abs(fmasks1[0].data()) ) > 4.0
  #
  xs.set_occupancies(value = 0)
  mp.ignore_zero_occupancy_atoms = False
  mp.use_asu_masks = False
  mm2 = masks.manager(miller_array   = f_calc,
                      xray_structure = xs,
                      mask_params    = mp)
  assert list(xs.scatterers().extract_occupancies()) == [0.0, 0.0]
  fmasks1 = mm1.shell_f_masks()
  assert len(fmasks1) == 1
  fmasks2 = mm2.shell_f_masks()
  assert len(fmasks2)==1
  assert flex.mean( flex.abs(fmasks1[0].data()) ) == \
         flex.mean( flex.abs(fmasks2[0].data()) )
  #
  mp.ignore_zero_occupancy_atoms = True
  mp.use_asu_masks = False
  mm3 = masks.manager(miller_array   = f_calc,
                      xray_structure = xs,
                      mask_params    = mp)
  assert list(xs.scatterers().extract_occupancies()) == [0.0, 0.0]
  fmasks3 = mm3.shell_f_masks()
  assert len(fmasks3)==1
  assert approx_equal(
    flex.abs(fmasks3[0].data()).min_max_mean().as_tuple(), (0.0, 0.0, 0.0))

def exercise_centrics(space_group_info, n_sites=10):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=(("O","N","C")*(n_sites//3+1))[:n_sites],
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
          ignore_zero_occupancy_atoms = False,
          solvent_radius=solvent_radius,
          shrink_truncation_radius=shrink_truncation_radius)
        f_mask = bulk_solvent_mask.structure_factors(miller_set=miller_set)
        centrics = f_mask.select_centric()
        if (centrics.indices().size() > 0):
          ideal = centrics.phase_transfer(centrics)
          assert flex.max(flex.abs(ideal.data() - centrics.data())) < 1.e-6

def structure_init2(site,site2,sg,cell):
  symmetry = crystal.symmetry(unit_cell=cell,
                              space_group_symbol=sg)
  structure = xray.structure(crystal_symmetry=symmetry)
  scatterer = xray.scatterer(
                 site = site,
                 u = 0.1,
                 occupancy = 1.0,
                 scattering_type = "C")
  scatterer2 = xray.scatterer(
                 site = site2,
                 u = 0.1,
                 occupancy = 1.0,
                 scattering_type = "C")
  structure.add_scatterer(scatterer)
  structure.add_scatterer(scatterer2)
  xyzf = flex.vec3_double()
  atmrad = flex.double()
  for scatterer in structure.scatterers():
    xyzf.append(list(scatterer.site))
    atmrad.append(van_der_waals_radii.vdw.table[scatterer.element_symbol()])
  assert xyzf.size() == atmrad.size()
  return structure, xyzf, atmrad

def tst2_run(angles, nspacing, af ):
  rad = van_der_waals_radii.vdw.table["C"]
  a = (2.0*rad)*2.0*af;
  pos1=(0.25,0.25,0.25)
  pos2=(0.25+0.9*2.0*rad/a,0.25,0.25)
  solvent_radius1 = rad*0.5
  solvent_radius2 = rad*3.1
  # shrink_truncation_radius = 1.3*solvent_radius2
  shrink_truncation_radius = 0.0
  cell = [a,a,a, angles[0], angles[1], angles[2] ]
  structure1, xyzf1, atmrad1 = structure_init2(pos1,pos2, "P1", cell)
  structure2, xyzf2, atmrad2 = structure_init2(pos1,pos2, "P1", cell)
  step = (1.0/nspacing) * a
  #
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell = structure1.unit_cell(),
    step = step)
  #
  m1 = masks.around_atoms(
    structure1.unit_cell(),
    structure1.space_group().order_z(),
    xyzf1,
    atmrad1,
    crystal_gridding.n_real(),
    solvent_radius1,
    shrink_truncation_radius, explicit_distance=False, debug=True )
  #
  m2 = masks.around_atoms(
    structure2.unit_cell(),
    structure2.space_group().order_z(),
    xyzf2,
    atmrad2,
    crystal_gridding.n_real(),
    solvent_radius2,
    shrink_truncation_radius, explicit_distance=False, debug=True )
  #
  m3 = masks.around_atoms(
    structure1.unit_cell(),
    structure1.space_group().order_z(),
    xyzf1,
    atmrad1,
    crystal_gridding.n_real(),
    solvent_radius1,
    shrink_truncation_radius, explicit_distance=True, debug=True)
  #
  m4 = masks.around_atoms(
    structure2.unit_cell(),
    structure2.space_group().order_z(),
    xyzf2,
    atmrad2,
    crystal_gridding.n_real(),
    solvent_radius2,
    shrink_truncation_radius, explicit_distance=True, debug=True)
  #
  assert m3.n_atom_points == m4.n_atom_points
  mx1 = max(abs(m1.n_atom_points),abs(m2.n_atom_points))
  mx2 = max(abs(m1.n_atom_points),abs(m3.n_atom_points))
  mx3 =  max(abs(m2.n_atom_points),abs(m3.n_atom_points))
  if mx1!=0:
    assert float(abs(m1.n_atom_points - m2.n_atom_points))/float(mx1) < 0.01
  if mx2!=0:
    assert float(abs(m1.n_atom_points - m3.n_atom_points))/float(mx2) < 0.01
  if mx3!=0:
    assert float(abs(m2.n_atom_points - m3.n_atom_points))/float(mx3) < 0.01

def exercise_4():
 if(1):
  tst2_run( [10.0,90.0,90.0], 139.0, 4.0 )
 if(1):
  tst2_run( [ 90.0,90.0,10.0], 139.0, 4.0 )
 if(1):
  tst2_run( [90.0,10.0,90.0], 139.0, 4.0 )
 if(1):
  tst2_run( [90.0, 90.0,170.0], 139.0, 5.0 )
 if(1):
  tst2_run( [20.0,30.0,40.0], 139.0, 4.0 )
 if(1):
  tst2_run( [90.0,170.0,90.0], 139.0, 4.0 )

def run_call_back(flags, space_group_info):
  exercise_centrics(space_group_info)

def run():
  exercise_1()
  exercise_2()
  exercise_3()
  exercise_4()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back,
    symbols_to_stdout=True, symbols_to_stderr=False)

if (__name__ == "__main__"):
  run()
