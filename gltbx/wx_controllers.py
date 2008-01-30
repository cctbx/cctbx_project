from __future__ import division
import wx
import wx.lib.colourselect
import math
from libtbx import adopt_init_args

class material(object):

  def __init__(self,
               model,
               front_colour_picker,
               back_colour_picker,
               swap_colour_btn,
               ambient_slider,
               diffuse_slider,
               specular_slider,
               specular_focus_slider,
               on_change=lambda:None):
    adopt_init_args(self, locals())

    for slider in (self.ambient_slider, self.diffuse_slider,
                   self.specular_slider):
      slider.SetRange(0,100)
    self.specular_focus_slider.SetRange(-10, 30)

    self.front_colour_picker.Bind(wx.lib.colourselect.EVT_COLOURSELECT,
                                  self.on_front_colour_changed)
    self.back_colour_picker.Bind(wx.lib.colourselect.EVT_COLOURSELECT,
                                 self.on_back_colour_changed)
    self.ambient_slider.Bind(wx.EVT_SCROLL,
                             self.on_ambient_changed)
    self.diffuse_slider.Bind(wx.EVT_SCROLL,
                             self.on_diffuse_changed)
    self.specular_slider.Bind(wx.EVT_SCROLL,
                              self.on_specular_changed)
    self.specular_focus_slider.Bind(wx.EVT_SCROLL,
                                    self.on_specular_focus_changed)
    self.swap_colour_btn.Bind(wx.EVT_BUTTON, self.on_front_back_colours_swap)

    self.on_model_changed()

  def on_front_back_colours_swap(self, e):
    self.model.front_colour, self.model.back_colour\
        = self.model.back_colour, self.model.front_colour
    self.on_model_changed()
    self.on_change()

  def on_front_colour_changed(self, e):
    c = self.front_colour_picker.GetColour()
    self.model.front_colour = (c.Red()/255, c.Green()/255, c.Blue()/255)
    self.on_change()

  def on_back_colour_changed(self, e):
    c = self.back_colour_picker.GetColour()
    self.model.back_colour = (c.Red()/255, c.Green()/255, c.Blue()/255)
    self.on_change()

  def on_ambient_changed(self, e):
    self.model.ambient = self.ambient_slider.GetValue()/100
    self.on_change()

  def on_diffuse_changed(self, e):
    self.model.diffuse = self.diffuse_slider.GetValue()/100
    self.on_change()

  def on_specular_changed(self, e):
    self.model.specular = self.specular_slider.GetValue()/100
    self.on_change()

  def on_specular_focus_changed(self, e):
    self.model.specular_focus = 10**(self.specular_focus_slider.GetValue()/10)
    self.on_change()

  def on_model_changed(self):
    r,g,b = [ x*255 for x in self.model.front_colour[0:3] ]
    self.front_colour_picker.SetColour(wx.Colour(r,g,b))
    r,g,b = [ x*255 for x in self.model.back_colour[0:3] ]
    self.back_colour_picker.SetColour(wx.Colour(r,g,b))
    self.ambient_slider.SetValue(self.model.ambient*100)
    self.diffuse_slider.SetValue(self.model.diffuse*100)
    self.specular_slider.SetValue(self.model.specular*100)
    self.specular_focus_slider.SetValue(
      int(math.log10(self.model.specular_focus*10)))
