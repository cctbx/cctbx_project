from cctbx.array_family import flex
import sys
from cctbx import xray
from cctbx import crystal
from mmtbx.max_lik import max_like_non_uniform
from cctbx.development import random_structure
from cctbx.development import debug_utils


#------------------------------------------------------------------------TEST-0
def test_less_one(space_group_info,
                  volume_per_atom = 100,
                  d_min = 1.8):

  #symmetry = crystal.symmetry(space_group_symbol="p212121")
  #space_group_info = symmetry.space_group_info()
  #n_sym = space_group_info.type().group().n_smx()
  for n_atoms in (1,10):
    if(n_atoms == 1):
       denom = 16.0
    else:
       denom = 4.0

    for element in ("C","O","N"):
        structure = random_structure.xray_structure(
               space_group_info=space_group_info,
               elements=[element]*n_atoms,
               volume_per_atom=volume_per_atom,
               random_u_iso=False)

        fc = structure.structure_factors(d_min          = d_min,
                                         anomalous_flag = False,
                                         algorithm      = "fft").f_calc()
        manager = max_like_non_uniform.ordered_solvent_distribution(
                                            structure = structure,
                                            fo        = fc,
                                            grid_step = fc.d_min()/denom)
        f_water_dist = manager.fcalc_from_distribution()
        ### check phase compatibility with the symmetry:
        centrics = f_water_dist.select_centric()
        if(centrics.indices().size() > 0):
           ideal = centrics.phase_transfer(centrics)
           assert flex.max(flex.abs(ideal.data() - centrics.data())) < 1.e-6
        ###
        #print "max = ", flex.max( flex.abs( f_water_dist.data() ) )
        #print "min = ", flex.min( flex.abs( f_water_dist.data() ) )
        #print "ave = ", flex.mean( flex.abs( f_water_dist.data() ) )
        assert flex.max( flex.abs( f_water_dist.data() ) ) < 1.0
#------------------------------------------------------------------------TEST-1
def test_grid_step(n_sites = 50,
                   volume_per_atom = 50,
                   d_min = 2.0):
  grid_step = (0.2,0.4,0.6,0.7,0.9,1.0)
  for step in grid_step:
    symmetry = crystal.symmetry(space_group_symbol="P1")
    structure = random_structure.xray_structure(space_group_info = symmetry.space_group_info(),
                                                elements=["C"]*n_sites,
                                                volume_per_atom=volume_per_atom,
                                                random_u_iso=False)
    fc = structure.structure_factors(d_min = d_min,
                                     anomalous_flag=False,
                                     algorithm="fft").f_calc()
    manager = max_like_non_uniform.ordered_solvent_distribution(
                                      structure = structure,
                                      fo = fc,
                                      grid_step = step)
    f_water_dist = manager.fcalc_from_distribution()
    ### check phase compatibility with the symmetry:
    centrics = f_water_dist.select_centric()
    if(centrics.indices().size() > 0):
       ideal = centrics.phase_transfer(centrics)
       assert flex.max(flex.abs(ideal.data() - centrics.data())) < 1.e-6
    ###
    #print "max = ", flex.max( flex.abs( f_water_dist.data() ) )
    #print "min = ", flex.min( flex.abs( f_water_dist.data() ) )
    #print "ave = ", flex.mean( flex.abs( f_water_dist.data() ) )
    assert flex.max( flex.abs( f_water_dist.data() ) ) < 1.0

#------------------------------------------------------------------------TEST-2
def test_r(space_group_info,
           step = 0.6,
           d_min = 4.0):
  r = (-1.05,-0.1,0.0,1.05,0.1)
  n_sym = space_group_info.type().group().n_smx()
  if(n_sym != 8):
     for x in r:
      for y in r:
       for z in r:
         symmetry = crystal.symmetry(unit_cell = space_group_info.any_compatible_unit_cell(volume= 2000),
                                     space_group_info = space_group_info)
         #symmetry = crystal.symmetry(unit_cell=(11., 12., 13., 75., 85., 95.),
         #                            space_group_symbol="P1")
         #symmetry = crystal.symmetry(unit_cell=(11., 12., 13., 90., 90., 90.),
         #                            space_group_symbol="p212121")
         structure = xray.structure(crystal_symmetry=symmetry)
         scatterer = xray.scatterer(
                          site = (x,y,z),
                          u = 0.1,
                          occupancy = 1.0,
                          scattering_type = "O")
         structure.add_scatterer(scatterer)
         fc = structure.structure_factors(d_min = d_min,
                                          anomalous_flag=False,
                                          algorithm="fft").f_calc()
         manager = max_like_non_uniform.ordered_solvent_distribution(
                                         structure = structure,
                                         fo = fc,
                                         grid_step = step)
         f_water_dist = manager.fcalc_from_distribution()
         ### check phase compatibility with the symmetry:
         centrics = f_water_dist.select_centric()
         if(centrics.indices().size() > 0):
            ideal = centrics.phase_transfer(centrics)
            assert flex.max(flex.abs(ideal.data() - centrics.data())) < 1.e-6
         ###
         #print "max = ", flex.max( flex.abs( f_water_dist.data() ) )
         #print "min = ", flex.min( flex.abs( f_water_dist.data() ) )
         #print "ave = ", flex.mean( flex.abs( f_water_dist.data() ) )
         assert flex.max( flex.abs( f_water_dist.data() ) ) < 1.0
#------------------------------------------------------------------------

def run_call_back(flags, space_group_info):
  test_less_one(space_group_info = space_group_info)

def run_call_back_1(flags, space_group_info):
  test_r(space_group_info = space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back_1,
    symbols_to_stdout=True, symbols_to_stderr=False)
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back,
    symbols_to_stdout=True, symbols_to_stderr=False)
  test_grid_step()
  print "OK"

if (__name__ == "__main__"):
  run()
