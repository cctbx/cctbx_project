from scitbx.iso_surface import tst_iso_surface
from crys3d import wx_map_viewer
from cctbx import uctbx


def exercise():
  uc = uctbx.unit_cell((1,1,1,60,120,90))
  #uc = uctbx.unit_cell((1,1,1,90,90,90))
  elliptic = tst_iso_surface.triangulation_test_case(
    func=tst_iso_surface.elliptic(),
    grid_size=(50, 40, 30),
    lazy_normals=False,
    descending_normals=True)
  wx_map_viewer.display(unit_cell=uc,
                        raw_map=elliptic.map,
                        iso_level=lambda map_stats: 1.3,
                        from_here=(-0.2, -0.3, -0.4),
                        to_there=(0.3, 0.4, 0.2),
                        #periodic=True,
                        wires=False,
                        title="Ellipsoid")

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
