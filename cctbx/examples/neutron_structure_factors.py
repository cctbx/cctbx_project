from __future__ import absolute_import, division, print_function
from cctbx import neutron
from cctbx import crystal

def run():
  quartz_structure = neutron.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(5.01,5.01,5.47,90,90,120),
        space_group_symbol="P6222")),
    scatterers=[
      neutron.scatterer(
        label="Si",
        site=(1/2.,1/2.,1/3.),
        u=0.2),
      neutron.scatterer(
        label="O",
        site=(0.197,-0.197,0.83333),
        u=(.01,0.01,0.01,0,0,0))])
  quartz_structure.show_summary().show_scatterers()
  f_calc = quartz_structure.structure_factors(d_min=1)
  f_calc.show_summary().show_array()

if (__name__ == "__main__"):
  run()
