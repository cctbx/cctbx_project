"""\
Updated example from IUCr Computing Commission Newsletter No. 1:
http://cci.lbl.gov/publications/download/iucrcompcomm_jan2003.pdf
"""

from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

def run():
  quartz_structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(5.01,5.01,5.47,90,90,120),
        space_group_symbol="P6222")),
    scatterers=flex.xray_scatterer([
      xray.scatterer(
        label="Si",
        site=(1/2.,1/2.,1/3.),
        u=0.2),
      xray.scatterer(
        label="O",
        site=(0.197,-0.197,0.83333),
        u=0)]))

  quartz_structure.show_summary().show_scatterers()

  from libtbx import easy_pickle
  easy_pickle.dump("beach", quartz_structure)

  from libtbx import easy_pickle
  quartz_structure = easy_pickle.load("beach")

  for scatterer in quartz_structure.scatterers():
    print "%s:" % scatterer.label, "%8.4f %8.4f %8.4f" % scatterer.site
    site_symmetry = quartz_structure.site_symmetry(scatterer.site)
    print "  point group type:", site_symmetry.point_group_type()
    print "  special position operator:", site_symmetry.special_op_simplified()

  for table in ["xray", "electron"]:
    print "Scattering type table:", table

    reg = quartz_structure.scattering_type_registry(table=table)
    reg.show_summary()

    f_calc = quartz_structure.structure_factors(d_min=2).f_calc()
    f_calc.show_summary().show_array()

    f_calc.d_spacings().show_array()

    low_resolution_only = f_calc.select(f_calc.d_spacings().data() > 2.5)
    low_resolution_only.show_array()

    print

if (__name__ == "__main__"):
  run()
