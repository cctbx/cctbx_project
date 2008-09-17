from scitbx.iso_surface import tst_iso_surface
from crys3d import wx_map_viewer
from cctbx import uctbx


def exercise():
  elliptic = tst_iso_surface.triangulation_test_case(
    func=tst_iso_surface.sinusoidal(),
    grid_size=(50, 40, 30),
    lazy_normals=False,
    descending_normals=True)
  class my_map_view(wx_map_viewer.map_view):
    def set_initial_iso_level(self, density_stats):
      self.iso_level = 0.15
  wx_map_viewer.display(unit_cell=uctbx.unit_cell((1,1,1,90,90,90)),
                        map_view_type=my_map_view,
                        raw_map=elliptic.map,
                        wires=False,
                        title="Ellipsoid")

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
