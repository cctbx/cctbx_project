from scitbx.iso_surface import tst_iso_surface
from crys3d import wx_map_viewer
from cctbx import uctbx


def exercise():
  uc = uctbx.unit_cell((1,1,1,60,120,90))
  #uc = uctbx.unit_cell((1,1,1,90,90,90))
  elliptic = tst_iso_surface.triangulation_test_case(
    func=tst_iso_surface.sinusoidal(),
    grid_size=(50, 40, 30),
    lazy_normals=False,
    descending_normals=True)
  wx_map_viewer.display(unit_cell=uc,
                        raw_map=elliptic.map,
                        iso_level=lambda map_stats: 0.15,
                        wires=False,
                        title="Ellipsoid")

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
