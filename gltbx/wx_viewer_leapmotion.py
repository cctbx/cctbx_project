
"""
Extensions for using the LeapMotion controller with the wxGLWindow class.
Should be safe to use even if no controller (and no supporting Python module)
is present.

The basic functionality (though not the actual code) is cribbed from Brad
Larson's BSD-licensed Molecules app for Mac.
Copyright (c) 2013, Sunset Lake Software
All rights reserved.
https://github.com/BradLarson/MoleculesMac
"""

from __future__ import absolute_import, division, print_function
import wx

leap_enabled = False
try :
  import Leap
except ImportError :
  Leap = None
else :
  leap_enabled = True
  #from Leap import CircleGesture, KeyTapGesture, ScreenTapGesture, SwipeGesture

LEAP_TRANSLATE_ID = wx.NewId()
LEAP_ROTATE_ID = wx.NewId()
LEAP_SCALE_ID = wx.NewId()

class LeapMotionEvent(wx.PyEvent):
  event_id = None
  def __init__(self, data, **kwds):
    self.data = data
    self.__dict__.update(kwds)
    wx.PyEvent.__init__(self)
    self.SetEventType(self.event_id)

class LeapTranslationEvent(LeapMotionEvent):
  event_id = LEAP_TRANSLATE_ID

class LeapRotationEvent(LeapMotionEvent):
  event_id = LEAP_ROTATE_ID

class LeapScaleEvent(LeapMotionEvent):
  event_id = LEAP_SCALE_ID

if (Leap is not None):
  class Listener(Leap.Listener):
    def __init__(self, viewer):
      self._prev_frame = None
      self.viewer = viewer
      Leap.Listener.__init__(self)

    def on_init(self, controller):
      #print "Initialized"
      self._prev_frame = None

    def on_connect(self, controller):
      #print "Connected"
      controller.enable_gesture(Leap.Gesture.TYPE_CIRCLE);
      controller.enable_gesture(Leap.Gesture.TYPE_KEY_TAP);
      controller.enable_gesture(Leap.Gesture.TYPE_SCREEN_TAP);
      controller.enable_gesture(Leap.Gesture.TYPE_SWIPE);

    def on_disconnect(self, controller):
      pass #print "Disconnected"

    def on_exit(self, controller):
      pass #print "Exited"

    def on_frame(self, controller):
      frame = controller.frame()
      prev_frame = self._prev_frame
      self._prev_frame = frame
      if (prev_frame is None):
        return
      open_hands = []
      for hand in frame.hands :
        if (len(hand.fingers) > 1):
          open_hands.append(hand)
      if (len(open_hands) < 1):
        self._prev_frame = None
      elif (len(open_hands) == 1):
        first_hand = open_hands[0]
        if (len(first_hand.fingers) > 2):
          translation = first_hand.translation(prev_frame)
          dx_abs = abs(translation.x)
          dy_abs = abs(translation.y)
          dz_abs = abs(translation.z)
          if (dx_abs > 40) or (dy_abs > 40) or (dz_abs > 40):
            self._prev_frame = None
            return
          elif (dz_abs > 1):
            self.scale(translation)
          elif (dx_abs > 1) or (dy_abs > 1):
            self.rotate(translation)
      else :
        translation = frame.translation(prev_frame)
        dx_abs = abs(translation.x)
        dy_abs = abs(translation.y)
        dz_abs = abs(translation.z)
        if (dx_abs > 40) or (dy_abs > 40) or (dz_abs > 40):
          self._prev_frame = None
          return
        if (dx_abs > 1) or (dy_abs > 1):
          self.translate(translation)

    def translate(self, vector):
      event = LeapTranslationEvent((vector[0], vector[1]))
      wx.PostEvent(self.viewer, event)

    def rotate(self, vector):
      event = LeapRotationEvent((vector[0], vector[1]))
      wx.PostEvent(self.viewer, event)

    def scale(self, vector):
      event = LeapScaleEvent(data=vector.z*0.02)
      wx.PostEvent(self.viewer, event)

else :
  class Listener(object):
    def __init__(self, viewer) : pass

class wxLeapMotionWindowMixin(object):
  def __init__(self):
    assert isinstance(self, wx.Window)
    self.Connect(-1, -1, LEAP_TRANSLATE_ID, self.OnLeapTranslate)
    self.Connect(-1, -1, LEAP_ROTATE_ID, self.OnLeapRotate)
    self.Connect(-1, -1, LEAP_SCALE_ID, self.OnLeapScale)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)
    self._leap_controller = None
    self._leap_listener = None

  def start_leapmotion_controller_listener(self):
    if (Leap is not None):
      self._leap_controller = Leap.Controller()
      self._leap_listener = Listener(self)
      self._leap_controller.add_listener(self._leap_listener)

  def stop_leapmotion_controller_listener(self):
    if (self._leap_controller is not None):
      self._leap_controller.remove_listener(self._leap_listener)

  def OnDestroy(self, event):
    self.stop_leapmotion_controller_listener()

  def OnLeapTranslate(self, event):
    pass

  def OnLeapRotate(self, event):
    pass

  def OnLeapScale(self, event):
    pass

#-----------------------------------------------------------------------
# CCTBX-SPECIFIC CODE
from gltbx.wx_viewer import wxGLWindow
import gltbx.util
from scitbx import matrix
class wxGLWindowLeapEnabled(wxGLWindow, wxLeapMotionWindowMixin):
  # TODO make these user-configurable
  _leap_rscale = 5
  _leap_tscale = 1/3.
  """
  Drop-in replacement for gltbx.wx_viewer.wxGLWindow.  The Listener object
  is created when self.initialize_modelview() is called.  Instead of altering
  the view directly, the OnLeap* event handlers set internal attributes that
  are processed as part of the idle loop, which seems to run much more
  smoothly.
  """
  def __init__(self, *args, **kwds):
    wxGLWindow.__init__(self, *args, **kwds)
    wxLeapMotionWindowMixin.__init__(self)
    self._leap_rotation = None
    self._leap_translation = None
    self._leap_scale = None

  def initialize_modelview(self):
    wxGLWindow.initialize_modelview(self)
    if (self._leap_controller is None):
      self.start_leapmotion_controller_listener()

  def OnLeapTranslate(self, event):
    self._leap_translation = event.data
    event.Skip()

  def OnLeapRotate(self, event):
    self._leap_rotation = event.data
    event.Skip()

  def OnLeapScale(self, event):
    self._leap_scale = event.data
    #self.OnScale(event.data)
    event.Skip()

  def OnIdle(self, event):
    if (self._leap_scale is not None):
      self.OnScale(self._leap_scale)
      self._leap_scale = None
    if (self._leap_rotation is not None):
      self.leap_rotate()
    if (self._leap_translation is not None):
      self.leap_translate()

  def leap_rotate(self):
    dx, dy = self._leap_rotation
    self.rotate_view(0, 0, dx*self._leap_rscale, -dy*self._leap_rscale,
      shift_down=False)
    self._leap_rotation = None

  def leap_translate(self):
    rc = self.rotation_center
    rc_eye = gltbx.util.object_as_eye_coordinates(rc)
    dx, dy = self._leap_translation
    gltbx.util.translate_object(1, 0, 0, -dx*self._leap_tscale,
      dy*self._leap_tscale)
    self._leap_translation = None
    self.rotation_center = tuple(
      gltbx.util.modelview_matrix_as_rt().inverse() * matrix.col(rc_eye))
    self.OnRedraw()
