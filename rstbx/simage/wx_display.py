import wx

class detector_surface(wx.Window):

  def __init__(O, parent, work_params):
    super(detector_surface, O).__init__(parent=parent, id=-1)
    O.work_params = work_params
    O.Bind(wx.EVT_SIZE, O.OnSize)
    O.Bind(wx.EVT_PAINT, O.OnPaint)
    O.Bind(wx.EVT_CHAR, O.OnChar)
    O.reset_state()

  def reset_state(O):
    O.prev_work_phil_ewp_none_str = None
    O.prev_work_phil_str = None
    O.pixels = None
    O.image = None
    O.image_2 = None
    O.spots = None
    O.predicted_spots = None
    O.image_toggle = False

  def recompute(O):
    w, h = O.GetSizeTuple()
    p = min(w, h)
    if (p == 0):
      O.reset_state()
      return False
    O.work_params.detector.pixels = (p, p)
    O.GetParent().wx_detector_pixels.SetLabel("%d x %d" % (p, p))
    ewp_save = O.work_params.ewald_proximity
    O.work_params.ewald_proximity = None
    work_phil_ewp_none_str = O.work_params.phil_master.format(
      O.work_params).as_str()
    O.work_params.ewald_proximity = ewp_save
    work_phil_str = O.work_params.phil_master.format(O.work_params).as_str()
    if (   O.prev_work_phil_str is None
        or O.prev_work_phil_str != work_phil_str):
      O.prev_work_phil_str = work_phil_str
      from rstbx.simage.create import compute_image
      O.pixels = compute_image(O.work_params)
      saturation = int(_.signal_max * _.saturation_level + 0.5)
      O.image = O.pixels.as_rgb_scale_string(
        rgb_scales_low=(1,1,1),
        rgb_scales_high=(0,0,0),
        saturation=saturation)
      if (O.work_params.wavelength_2 is not None):
        pixels_2 = compute_image(O.work_params, use_wavelength_2=True)
        O.image_2 = pixels_2.as_rgb_scale_string(
          rgb_scales_low=(1,1,1),
          rgb_scales_high=(1,0,0),
          saturation=saturation)
      if (   O.prev_work_phil_ewp_none_str is None
          or O.prev_work_phil_ewp_none_str != work_phil_ewp_none_str):
        O.prev_work_phil_ewp_none_str = work_phil_ewp_none_str
        O.spots = None
        O.predicted_spots = None
    return True

  def run_spotfinder(O):
    if (O.pixels is None):
      return
    if (O.spots is not None):
      return
    dpx,dpy = O.work_params.detector.pixels
    if (dpx < 100 or dpy < 100):
      return
    from rstbx.simage import run_spotfinder
    O.spots = run_spotfinder.process(
      work_params=O.work_params,
      pixels=O.pixels,
      show_spots=False)
    O.Refresh()

  def run_labelit_index(O, use_original_uc_cr=False):
    if (O.predicted_spots is not None):
      return
    O.run_spotfinder()
    if (O.spots is None):
      return
    if (O.spots.size() < 10):
      print "Insufficient number of spotfinder spots."
      print
      return
    else:
      if (use_original_uc_cr):
        print
        print "Using original unit cell and crystal rotation" \
          " for spot prediction."
        print
        uc = O.work_params.unit_cell
        cr = O.work_params.crystal_rotation_matrix
      else:
        from rstbx.simage.run_labelit_index import process
        try:
          ai = process(work_params=O.work_params, spots=O.spots)
        except Exception, e:
          print "Indexing exception:", e
          print
          return
        print "Spots indexed: %d of %d" % (
          ai.hklobserved().size(),
          O.spots.size())
        co = ai.getOrientation()
        uc = co.unit_cell()
        cr = co.crystal_rotation_matrix()
        print "labelit index unit cell:", uc
        from rstbx.simage import refine_uc_cr
        refined = refine_uc_cr.refine(
          work_params=O.work_params,
          spots=O.spots,
          good_i_seqs=ai.get_observed_spot_i_seqs_good_only(),
          miller_indices=ai.hklobserved(),
          unit_cell=uc,
          crystal_rotation=cr)
        uc = refined.unit_cell
        cr = refined.crystal_rotation
        print
      import cctbx.crystal
      crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc,
        space_group_symbol="P1")
      d_min = O.work_params.d_min
      if (d_min is None):
        d_min = O.work_params.wavelength
      import cctbx.miller
      miller_set = cctbx.miller.build_set(
        crystal_symmetry=crystal_symmetry,
        d_min=d_min,
        anomalous_flag=True)
      from rstbx.simage import image_simple
      O.predicted_spots = image_simple(store_spots=True).compute(
        unit_cell=miller_set.unit_cell(),
        miller_indices=miller_set.indices(),
        spot_intensity_factors=None,
        crystal_rotation_matrix=cr,
        ewald_radius=1/O.work_params.wavelength,
        ewald_proximity=O.work_params.ewald_proximity,
        signal_max=O.work_params.signal_max,
        detector_distance=O.work_params.detector.distance,
        detector_size=O.work_params.detector.size,
        detector_pixels=O.work_params.detector.pixels,
        point_spread=O.work_params.point_spread,
        gaussian_falloff_scale=O.work_params.gaussian_falloff_scale).spots
      print "Number of predicted spots:", O.predicted_spots.size()
      print
    O.Refresh()

  def draw_image(O):
    w, h = O.work_params.detector.pixels
    assert w == h
    p = w
    assert p != 0
    wx_image = wx.EmptyImage(p, p)
    if (O.image is not None):
      if (not O.image_toggle): im = O.image
      else:                    im = O.image_2
      wx_image.SetData(im)
    wx_bitmap = wx_image.ConvertToBitmap()
    dc = wx.PaintDC(win=O)
    dc = wx.BufferedDC(dc)
    dc.Clear()
    w, h = O.GetSizeTuple()
    dc.DrawBitmap(bmp=wx_bitmap, x=w-p, y=0, useMask=False)
    dc.SetPen(wx.Pen("GREY", 0))
    dc.SetBrush(wx.Brush("GREY", wx.CROSSDIAG_HATCH))
    if (w > p):
      dc.DrawRectangle(0, 0, w-(p+1), h)
    elif (h > p):
      dc.DrawRectangle(0, p+1, w, h-(p+1))
    if (O.spots is not None):
      dc.SetPen(wx.Pen("RED", 1))
      dc.SetBrush(wx.Brush("WHITE", wx.TRANSPARENT))
      for spot in O.spots:
        x,y = spot.ctr_mass_x()+0.5, spot.ctr_mass_y()+0.5
        dc.DrawCircle(x=w-p+y, y=x, radius=5)
    if (O.predicted_spots is not None):
      dc.SetPen(wx.Pen("BLUE", 1))
      dc.SetBrush(wx.Brush("WHITE", wx.TRANSPARENT))
      for spot in O.predicted_spots:
        x,y,_ = spot
        dc.DrawCircle(x=w-p+y, y=x, radius=5)

  def OnSize(O, event):
    if (O.recompute()):
      O.Refresh()

  def OnPaint(O, event):
    if (O.image is None and not O.recompute()):
      return
    O.draw_image()

  def OnChar(O, event):
    key = event.GetKeyCode()
    if (key == ord("w")):
      if (O.image_2 is not None):
        O.image_toggle = not O.image_toggle
        O.Refresh()
    elif (key == ord("s")):
      O.run_spotfinder()
    elif (key == ord("i")):
      O.predicted_spots = None
      O.run_labelit_index()
    elif (key == ord("I")):
      O.predicted_spots = None
      O.run_labelit_index(use_original_uc_cr=True)
    else:
      print "No action for this key stroke."

class main_panel(wx.Panel):

  def __init__(O, parent, work_params):
    super(main_panel, O).__init__(parent=parent, id=-1)
    O.work_params = work_params

    O.variable_name_by_wx_id = {}
    O.variable_values_by_name = {}
    v_sizer = wx.BoxSizer(orient=wx.VERTICAL)
    def add_slider(variable_name, val_min_max):
      h_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
      ctrl_id = wx.NewId()
      slider = wx.Slider(
        parent=O,
        id=ctrl_id,
        value=val_min_max[0],
        minValue=val_min_max[1],
        maxValue=val_min_max[2],
        size=(150, -1),
        style=wx.SL_HORIZONTAL|wx.SL_AUTOTICKS|wx.SL_LABELS)
      O.variable_name_by_wx_id[ctrl_id] = variable_name
      O.variable_values_by_name[variable_name] = val_min_max[0]
      slider.Bind(event=wx.EVT_SCROLL, handler=O.OnSliderScroll)
      h_sizer.Add(item=slider)
      label = wx.StaticText(parent=O, id=-1, label=variable_name)
      h_sizer.Add(item=label)
      v_sizer.Add(item=h_sizer)
      v_sizer.AddSpacer(3)
    ucp = O.work_params.unit_cell.parameters()
    for variable_name,value in zip(["a", "b", "c"], ucp[:3]):
      add_slider(variable_name, (
        value,
        min(10, round(value-0.5, 0)),
        max(100, round(value+5, -1))))
    for variable_name,value in zip(["alpha", "beta", "gamma"], ucp[3:]):
      add_slider(variable_name, (
        value,
        min(60, round(value-5, -1)),
        max(105, round(value+5, -1))))
    from libtbx.math_utils import normalize_angle
    xyz = O.work_params.euler_angles_xyz
    for variable_name,value in zip(["rot x", "rot y", "rot z"], xyz):
      add_slider(variable_name, (
        normalize_angle(value, deg=True, zero_centered=True),
        -180,
        180))

    def add_fs(min_val, max_val, increment, digits, label, value):
      import wx.lib.agw.floatspin as FS
      ctrl_id = wx.NewId()
      fs = FS.FloatSpin(
        parent=O,
        id=ctrl_id,
        min_val=min_val,
        max_val=max_val,
        increment=increment,
        value=value,
        agwStyle=FS.FS_LEFT)
      O.variable_name_by_wx_id[ctrl_id] = label
      fs.SetFormat("%f")
      fs.SetDigits(digits)
      fs.Bind(event=FS.EVT_FLOATSPIN, handler=O.OnFloatSpin),
      h_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
      h_sizer.Add(item=fs)
      label = wx.StaticText(parent=O, id=-1, label=label)
      h_sizer.Add(item=label, flag=wx.LEFT, border=5)
      v_sizer.Add(item=h_sizer)
      v_sizer.AddSpacer(3)

    add_fs(
      min_val=0.1,
      max_val=10,
      increment=0.1,
      digits=6,
      label="Wavelength",
      value=O.work_params.wavelength)

    add_fs(
      min_val=0.1,
      max_val=10,
      increment=0.1,
      digits=6,
      label="Wavelength 2",
      value=O.work_params.wavelength_2)

    if (O.work_params.d_min is None):
      O.work_params.d_min = O.work_params.wavelength
    add_fs(
      min_val=0.1,
      max_val=10,
      increment=0.1,
      digits=6,
      label="d-min",
      value=O.work_params.d_min)

    add_fs(
      min_val=-1,
      max_val=1,
      increment=0.01,
      digits=6,
      label="Ewald proximity",
      value=O.work_params.ewald_proximity)

    add_fs(
      min_val=1,
      max_val=500,
      increment=50,
      digits=2,
      label="Detector distance",
      value=O.work_params.detector.distance)

    def add_sp(label, value):
      sp = wx.SpinCtrl(parent=O, id=-1, min=1, max=100, initial=value)
      sp.Bind(event=wx.EVT_SPINCTRL, handler=O.OnSpin),
      h_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
      h_sizer.Add(item=sp)
      label = wx.StaticText(parent=O, id=-1, label=label)
      h_sizer.Add(item=label, flag=wx.LEFT, border=5)
      v_sizer.Add(item=h_sizer)
      v_sizer.AddSpacer(3)

    add_sp(
      label="Point spread",
      value=O.work_params.point_spread)

    add_fs(
      min_val=0,
      max_val=100,
      increment=1,
      digits=1,
      label="Gaussian falloff scale",
      value=O.work_params.gaussian_falloff_scale)

    def add_detector_pixels():
      h_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
      label = wx.StaticText(parent=O, id=-1, label="Detector pixels:")
      h_sizer.Add(item=label)
      dp = wx.StaticText(parent=O, id=-1, label="None")
      f = dp.GetFont()
      f.SetWeight(wx.BOLD)
      dp.SetFont(f)
      h_sizer.Add(item=dp, flag=wx.LEFT, border=5)
      v_sizer.AddSpacer(3)
      v_sizer.Add(item=h_sizer)
      return dp

    O.wx_detector_pixels = add_detector_pixels()

    O.detector_surface = detector_surface(parent=O, work_params=O.work_params)

    topsizer = wx.BoxSizer(orient=wx.HORIZONTAL)
    topsizer.Add(item=O.detector_surface, proportion=1, flag=wx.EXPAND)
    topsizer.Add(item=v_sizer, flag=wx.ALL, border=10)
    O.SetSizer(topsizer)
    topsizer.Layout()

  def reset_work_params(O):
    def getvar():
      return O.variable_values_by_name[variable_name]
    uc_params = []
    for variable_name in ["a", "b", "c", "alpha", "beta", "gamma"]:
      uc_params.append(getvar())
    from cctbx import uctbx
    try:
      uc = uctbx.unit_cell(uc_params)
    except RuntimeError:
      pass # simply keep previous
    else:
      O.work_params.unit_cell = uc
    xyz = []
    for variable_name in ["rot x", "rot y", "rot z"]:
      xyz.append(getvar())
    O.work_params.euler_angles_xyz = xyz

  def OnSliderScroll(O, event):
    variable_name = O.variable_name_by_wx_id[event.GetId()]
    O.variable_values_by_name[variable_name] = event.GetPosition()
    O.reset_work_params()
    O.detector_surface.recompute()
    O.detector_surface.Refresh()

  def OnFloatSpin(O, event):
    val = event.GetEventObject().GetValue()
    label = O.variable_name_by_wx_id[event.GetId()]
    if (label == "Wavelength"):
      O.work_params.wavelength = val
    elif (label == "Wavelength 2"):
      O.work_params.wavelength_2 = val
    elif (label == "d-min"):
      O.work_params.d_min = val
    elif (label == "Ewald proximity"):
      O.work_params.ewald_proximity = val
    elif (label == "Detector distance"):
      O.work_params.detector.distance = val
    elif (label == "Gaussian falloff scale"):
      O.work_params.gaussian_falloff_scale = val
    else:
      raise RuntimeError("Unknown label: %s" % label)
    O.detector_surface.recompute()
    O.detector_surface.Refresh()

  def OnSpin(O, event):
    val = event.GetEventObject().GetValue()
    O.work_params.point_spread = val
    O.detector_surface.recompute()
    O.detector_surface.Refresh()

def run(args):
  from rstbx.simage import run_spotfinder
  work_params = run_spotfinder.process_args(
    args=args, extra_phil_str="""\
saturation_level = 1.0
  .type = float
""")
  if (work_params.wavelength_2 is None):
    work_params.wavelength_2 = work_params.wavelength
  app = wx.App()
  frame = wx.Frame(
    parent=None,
    id=-1,
    title="wx_simage",
    pos=wx.DefaultPosition,
    size=wx.Size(800, 600))
  main_panel(parent=frame, work_params=work_params)
  frame.Show()
  app.MainLoop()
