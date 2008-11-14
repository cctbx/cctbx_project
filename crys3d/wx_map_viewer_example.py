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
  a = wx_map_viewer.App(unit_cell=uc,
                        raw_map=elliptic.map,
                        iso_level=lambda map_stats: 1.4,
                        from_here=(-0.3, 0.2, 0.4),
                        to_there=(0.7, 0.8, 1.6),
                        periodic=True,
                        wires=False,
                        title="Ellipsoid")
  #a.view_objects.orthographic = True
  a.MainLoop()

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
