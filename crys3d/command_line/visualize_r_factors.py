
# LIBTBX_SET_DISPATCHER_NAME cctbx.visualize_r_factors
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function
from crys3d.hklview.frames import *
from cctbx.miller.display import master_phil
from wxtbx import icons
from libtbx.utils import Sorry
import wx
import sys

class RfactorFrame (HKLViewFrame) :
  def add_view_specific_functions (self) :
    pass

  def load_reflections_file (self, file_name, **kwds) :
    if (isinstance(file_name, unicode)) :
      file_name = str(file_name)
    if (file_name != "") :
      from iotbx.reflection_file_reader import any_reflection_file
      from cctbx import miller
      from scitbx.array_family import flex
      try :
        hkl_file = any_reflection_file(file_name)
      except Exception as e :
        raise Sorry(str(e))
      arrays = hkl_file.as_miller_arrays(merge_equivalents=True)
      f_obs = f_model = None
      for array in arrays :
        labels = array.info().label_string()
        if labels.startswith("F-obs-filtered") :
          f_obs = array
        elif labels.startswith("F-model") :
          f_model = array
      if (f_obs is None) or (f_model is None) :
        raise Sorry("This does not appear to be a phenix.refine output "+
          "file.  The MTZ file should contain data arrays for the filtered "+
          "amplitudes (F-obs) and F-model.")
      f_delta = f_obs.customized_copy(sigmas=None,
        data=flex.abs(f_obs.data()-abs(f_model).data())).set_info(
          miller.array_info(labels=["abs(F_obs - F_model)"]))
      self.set_miller_array(f_delta)

def run (args) :
  import iotbx.phil
  pcl = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    reflection_file_def="data",
    usage_string="""\
cctbx.visualize_r_factors refine_001.mtz

Display the absolute difference between scaled and filtered F-obs and F-model,
which are used to calculate R-factors.
""")
  settings = pcl.work.extract()
  a = wx.App(0)
  app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_3d.GetBitmap())
  if (wx.VERSION >= (2,9)) :
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
  else :
    tb_icon = wx.TaskBarIcon()
  tb_icon.SetIcon(app_icon, "PHENIX data viewer")
  a.hklview_settings = settings
  f = RfactorFrame(None, -1, "F-model versus F-obs", size=(1024,768))
  f.Show()
  if (settings.data is not None) :
    f.load_reflections_file(settings.data)
  else :
    f.OnLoadFile(None)
  a.SetTopWindow(f)
  a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])
