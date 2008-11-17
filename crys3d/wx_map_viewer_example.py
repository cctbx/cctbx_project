from scitbx.iso_surface import tst_iso_surface
from crys3d import wx_map_viewer
from cctbx import uctbx


def exercise():
  uc = uctbx.unit_cell((1,1,1,60,120,90))
  #uc = uctbx.unit_cell((1,1,1,90,90,90))
  case = tst_iso_surface.triangulation_test_case(
    #func=tst_iso_surface.elliptic(),
    func=tst_iso_surface.periodic(),
    grid_size=(50, 40, 30),
    periodic=True,
    lazy_normals=False,
    descending_normals=True)
  a = wx_map_viewer.App(unit_cell=uc,
                        raw_map=case.map,
                        #iso_level=lambda map_stats: 3,
                        iso_level=lambda map_stats: 0.3,
                        #from_here=None, to_there=None,
                        from_here=(-0.5, -0.5, -0.5), to_there=(1.5, 1.5, 1.5),
                        wires=False,
                        title="Ellipsoid")
  #a.view_objects.orthographic = True
  a.MainLoop()

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
