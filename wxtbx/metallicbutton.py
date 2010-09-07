
# derived from wx.lib.platebtn (see copyright below).

###############################################################################
# Name: platebtn.py                                                           #
# Purpose: PlateButton is a flat label button with support for bitmaps and    #
#          drop menu.                                                         #
# Author: Cody Precord <cprecord@editra.org>                                  #
# Copyright: (c) 2007 Cody Precord <staff@editra.org>                         #
# Licence: wxWindows Licence                                                  #
###############################################################################

__all__ = ["MetallicButton", "AdjustAlpha", "AdjustColour",
           "GetHighlightColour",
           "GRADIENT_NORMAL", "GRADIENT_PRESSED", "GRADIENT_HIGHLIGHT",
           "MB_STYLE_DEFAULT", "GB_STYLE_BOLD_LABEL", "GB_STYLE_DROPARROW"]


import wx
import wx.lib.newevent
import wx.lib.wordwrap
import wx.lib.imageutils
from wx.lib.colourutils import *

# Used on OSX to get access to carbon api constants
if wx.Platform == '__WXMAC__':
    import Carbon.Appearance

GRADIENT_NORMAL = 0
GRADIENT_PRESSED = 1
GRADIENT_HIGHLIGHT = 2

MB_STYLE_DEFAULT = 1
MB_STYLE_BOLD_LABEL = 2
MB_STYLE_DROPARROW = 4

class MetallicButton (wx.PyControl) :
  def __init__ (self, parent, id_=wx.ID_ANY, label='', label2='', bmp=None,
                pos=wx.DefaultPosition, size=wx.DefaultSize,
                style=MB_STYLE_DEFAULT, name=wx.ButtonNameStr,
                start_color=(218,218,218), gradient_percent=15.0,
                highlight_color=(230,230,230), label_size=13,
                caption_size=11, button_margin=2, disable_after_click=0):
    wx.PyControl.__init__(self, parent, id_, pos, size,
      wx.BORDER_SIMPLE, name=name)
    self.InheritAttributes()
    self._bmp = dict(enable=bmp)
    self._margin = button_margin
    if bmp is not None :
      img = bmp.ConvertToImage()
      #img = img.ConvertToGreyscale(.795, .073, .026)
      wx.lib.imageutils.grayOut(img)
      self._bmp['disable'] = img.ConvertToBitmap()
    else :
      self._bmp['disable'] = None
    self._label2_font = self.GetFont()
    self._label2_font.SetPointSize(caption_size)
    self._label2_font.SetStyle(wx.FONTSTYLE_ITALIC)
    if style & MB_STYLE_BOLD_LABEL :
      font_size = label_size
      self._label_font = self.GetFont()
      self._label_font.SetPointSize(label_size)
      self._label_font.SetWeight(wx.FONTWEIGHT_BOLD)
      self.SetFont(self._label_font)
    #self._label2_font = wx.Font(caption_size, wx.SWISS, wx.ITALIC, wx.NORMAL)

    self._menu = None
    self.SetLabel(label)
    self._label2 = label2
    self._style = style
    self._size = tuple(size)
    self._state = dict(pre=GRADIENT_NORMAL, cur=GRADIENT_NORMAL)
    self._color = self.__InitColors(start_color, highlight_color,
      gradient_percent)
    self._caption_lines = None
    self._disable_after_click = disable_after_click

    # Setup Initial Size
    self.SetInitialSize()

    # Event Handlers
    self.Bind(wx.EVT_PAINT, lambda evt: self.__DrawButton())
    self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnErase)
    self.Bind(wx.EVT_SET_FOCUS, self.OnFocus)
    self.Bind(wx.EVT_KILL_FOCUS, self.OnKillFocus)

    # Mouse Events
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_LEFT_DCLICK, lambda evt: self.ToggleState())
    self.Bind(wx.EVT_ENTER_WINDOW,
              lambda evt: self.SetState(GRADIENT_HIGHLIGHT))
    self.Bind(wx.EVT_LEAVE_WINDOW,
              lambda evt: self.SetState(GRADIENT_NORMAL)) #wx.CallLater(10, self.SetState, GRADIENT_NORMAL))

    # Other events
    self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)
    self.Bind(wx.EVT_CONTEXT_MENU, lambda evt: self.ShowMenu())

  def __DrawBitmap(self, gc):
    """Draw the bitmap if one has been set
    @param gc: GCDC to draw with
    @return: x cordinate to draw text at

    """
    if self.IsEnabled():
      bmp = self._bmp['enable']
    else:
      bmp = self._bmp['disable']

    if bmp is not None and bmp.IsOk():
      bw, bh = bmp.GetSize()
      cw, ch = self.GetSize()
      if ch > (bh + 4) : # (self._margin * 2)):
        ypos = ((ch - bh) / 2) - (self._margin / 2)
        xpos = self._margin
      else :
        ypos = 0
        xpos = 0
      gc.DrawBitmap(bmp, xpos, ypos, bmp.GetMask() != None)
      return bw + 6
    else:
      return 6

  def __DrawDropArrow(self, gc, xpos, ypos):
    """Draw a drop arrow if needed and restore pen/brush after finished
    @param gc: GCDC to draw with
    @param xpos: x cord to start at
    @param ypos: y cord to start at

    """
    if self._menu is not None or self._style & MB_STYLE_DROPARROW:
      # Positioning needs a little help on Windows
      if wx.Platform == '__WXMSW__':
          xpos -= 2
      tripoints = [(xpos, ypos), (xpos + 6, ypos), (xpos + 3, ypos + 5)]
      brush_b = gc.GetBrush()
      pen_b = gc.GetPen()
      gc.SetPen(wx.TRANSPARENT_PEN)
      gc.SetBrush(wx.Brush(gc.GetTextForeground()))
      gc.DrawPolygon(tripoints)
      gc.SetBrush(brush_b)
      gc.SetPen(pen_b)
    else:
      pass

  def __DrawHighlight(self, gc, width, height):
    """Draw the main highlight/pressed state
    @param gc: GCDC to draw with
    @param width: width of highlight
    @param height: height of highlight

    """
    if self._state['cur'] == GRADIENT_PRESSED:
      color = self._color['press_start']
      end_color = self._color['press_end']
    else:
      color = self._color['hlight_start']
      end_color = self._color['hlight_end']

    rad = 0

    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    rgc = gc.GetGraphicsContext()
    brush = rgc.CreateLinearGradientBrush(0, 1, 0, height, color, end_color)
    rgc.SetBrush(brush)
    gc.DrawRectangle(0, 0, width, height)

  def __DrawCaption (self, gc, xpos, ypos) :
    if self._label2 != '' :
      gc.SetFont(self._label2_font)
      min_w, min_h = self._size
      if min_w == -1 :
        min_w = 120
      txt_w = min_w - xpos - 10
      if self._caption_lines is not None :
        lines = self._caption_lines
      else :
        lines = wx.lib.wordwrap.wordwrap(self._label2,
          width=min_w,
          dc=gc)
      offset = 0
      for line in lines.splitlines() :
        line_w, line_h = gc.GetTextExtent(line)
        gc.DrawText(line.rstrip(), xpos, ypos + offset)
        offset += line_h + 2

  def __PostEvent(self):
    """Post a button event to parent of this control"""
    bevt = wx.CommandEvent(wx.wxEVT_COMMAND_BUTTON_CLICKED, self.GetId())
    bevt.SetEventObject(self)
    bevt.SetString(self.GetLabel())
    wx.PostEvent(self.GetParent(), bevt)

  def __DrawButton(self):
    """Draw the button"""
    dc = wx.AutoBufferedPaintDCFactory(self)
    gc = wx.GCDC(dc)

    # Setup
    dc.SetBrush(wx.WHITE_BRUSH)
    gc.SetBrush(wx.WHITE_BRUSH)
    gc.SetFont(self.GetFont())
    #gc.SetBackgroundMode(wx.TRANSPARENT)

    # Calc Object Positions
    width, height = self.GetSize()
    tw, th = gc.GetTextExtent(self.GetLabel())
    if self._label2 != '' :
      txt_y = 4 #th + 4 #height - th - 4
      txt2_y = th + 8
    else :
      txt_y = max((height - th) / 2, 1)
      txt2_y = None
    #print height, th, txt_y, txt2_y

    if self._state['cur'] == GRADIENT_HIGHLIGHT:
      gc.SetTextForeground(self._color['htxt'])
      gc.SetPen(wx.TRANSPARENT_PEN)
      self.__DrawHighlight(gc, width, height)

    elif self._state['cur'] == GRADIENT_PRESSED:
      gc.SetTextForeground(self._color['htxt'])
      if wx.Platform == '__WXMAC__':
        brush = wx.Brush(wx.BLACK)
        brush.MacSetTheme(Carbon.Appearance.kThemeBrushFocusHighlight)
        pen = wx.Pen(brush.GetColour(), 1, wx.SOLID)
      else:
        pen = wx.Pen(AdjustColour(self._color['press_start'], -80, 220), 1)
      gc.SetPen(pen)

      self.__DrawHighlight(gc, width, height)
      txt_x = self.__DrawBitmap(gc)
      gc.DrawText(self.GetLabel(), txt_x + 2, txt_y)
      self.__DrawCaption(gc, txt_x + 2, txt2_y)
      self.__DrawDropArrow(gc, width - 10, (height / 2) - 2)

    else:
      rgc = gc.GetGraphicsContext()
      gc.SetPen(wx.TRANSPARENT_PEN)
      color =  wx.Colour(218,218,218)
      brush = rgc.CreateLinearGradientBrush(0, 1, 0, height,
        self._color['gradient_start'], self._color['gradient_end'])
      rgc.SetBrush(brush)
      gc.DrawRectangle(0, 0, width, height)
      if self.IsEnabled():
        gc.SetTextForeground(self.GetForegroundColour())
      else:
        txt_c = wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT)
        gc.SetTextForeground(txt_c)

    # Draw bitmap and text
    if self._state['cur'] != GRADIENT_PRESSED:
      txt_x = self.__DrawBitmap(gc)
      gc.DrawText(self.GetLabel(), txt_x + 2, txt_y)
      self.__DrawCaption(gc, txt_x + 2, txt2_y)
      #self.__DrawDropArrow(gc, txt_x + tw + 6, (height / 2) - 2)

  def __InitColors(self, start_color, highlight_color, gradient_percent):
    """Initialize the default colors"""
    start_color = wx.Colour(*start_color)
    start_hcolor = wx.Colour(*highlight_color) #GetHighlightColour()
    start_pcolor = AdjustColour(start_hcolor, -12)
    if gradient_percent != 0 :
      end_color = AdjustColour(start_color, gradient_percent)
      end_hcolor = AdjustColour(start_hcolor, gradient_percent)
      end_pcolor = AdjustColour(start_pcolor, gradient_percent)
    else :
      end_color = start_color
      end_hcolor = start_hcolor
      end_pcolor = start_pcolor
    colors = dict(default=True,
                  gradient_start=start_color,
                  gradient_end=end_color,
                  hlight_start=start_hcolor,
                  hlight_end=end_hcolor,
                  press_start=start_pcolor,
                  press_end=end_pcolor,
                  htxt=wx.Colour(0,0,0))
    # BestLabelColour(self.GetForegroundColour()))
    return colors

  #---- End Private Member Function ----#

  #---- Public Member Functions ----#
  def AcceptsFocus(self):
    """Can this window have the focus?"""
    return self.IsEnabled()

  @property
  def BitmapDisabled(self):
    """Property for accessing the bitmap for the disabled state"""
    return self._bmp['disable']

  @property
  def BitmapLabel(self):
    """Property for accessing the default bitmap"""
    return self._bmp['enable']

  # Aliases
  BitmapFocus = BitmapLabel
  BitmapHover = BitmapLabel
  BitmapSelected = BitmapLabel

  def Disable(self):
    """Disable the control"""
    wx.PyControl.Disable(self)
    self.Refresh()

  def DoGetBestSize(self):
    """Calculate the best size of the button
    @return: wx.Size

    """
    width = 8
    height = 10
    label_width = 0
    label_height = 0
    caption_width = 0
    caption_height = 0
    if self._bmp['enable'] is not None:
      bsize = self._bmp['enable'].GetSize()
      width += (bsize[0] + 10)
      height = bsize[1] + (self._margin * 2)
    else:
      width += 10

    if self.GetLabel():
      lsize = self.GetTextExtent(self.GetLabel())
      label_width = lsize[0]
      label_height = lsize[1]

    if self._label2 != '' :
      dc = wx.ClientDC(self)
      dc.SetFont(self._label2_font)
      min_w, min_h = self._size
      if min_w == -1 :
        min_w = 120
      txt_w = min_w - width - 10
      if wx.Platform == '__WXGTK__' :
        txt_w += 200
      lines = wx.lib.wordwrap.wordwrap(self._label2,
        width=txt_w,
        dc=dc)
      self._caption_lines = lines
      offset = 0
      if wx.Platform == "__WXMAC__" :
        buffer = 4
      else :
        buffer = 0
      for line in lines.splitlines() :
        line_w, line_h = dc.GetTextExtent(line)
        if line_w > caption_width :
          caption_width = line_w
        caption_height += line_h + buffer
    width += max(caption_width, label_width) + 4
    height = max(caption_height + label_height + 10, height)

    if self._menu is not None or self._style & MB_STYLE_DROPARROW :
       width += 12

    if width < self._size[0] :
      width = self._size[0]
    best = wx.Size(width, height)
    self.CacheBestSize(best)
    return best

  def Enable(self, enable=True):
    """Enable/Disable the control"""
    wx.PyControl.Enable(self, enable)
    self.Refresh()

  def GetBackgroundBrush(self, dc):
    """Get the brush for drawing the background of the button
    @return: wx.Brush
    @note: used internally when on gtk

    """
    if wx.Platform == '__WXMAC__' : #or self._style & PB_STYLE_NOBG:
      return wx.TRANSPARENT_BRUSH

    bkgrd = self.GetBackgroundColour()
    brush = wx.Brush(bkgrd, wx.SOLID)
    my_attr = self.GetDefaultAttributes()
    p_attr = self.GetParent().GetDefaultAttributes()
    my_def = bkgrd == my_attr.colBg
    p_def = self.GetParent().GetBackgroundColour() == p_attr.colBg
    if my_def and not p_def:
      bkgrd = self.GetParent().GetBackgroundColour()
      brush = wx.Brush(bkgrd, wx.SOLID)
    return brush

  def GetBitmapDisabled(self):
    """Get the bitmap of the disable state
    @return: wx.Bitmap or None

    """
    return self._bmp['disable']

  def GetBitmapLabel(self):
    """Get the label bitmap
    @return: wx.Bitmap or None

    """
    return self._bmp['enable']

  # GetBitmap Aliases for BitmapButton api
  GetBitmapFocus = GetBitmapLabel
  GetBitmapHover = GetBitmapLabel

  # Alias for GetLabel
  GetLabelText = wx.PyControl.GetLabel

  def GetMenu(self):
    """Return the menu associated with this button or None if no
    menu is associated with it.

    """
    return getattr(self, '_menu', None)

  def HasTransparentBackground(self):
    """Override setting of background fill"""
    return True

  @property
  def LabelText(self):
    """Property for getting the label of the button"""
    return self.GetLabel()

  #---- Event Handlers ----#

  def OnErase(self, evt):
    """Trap the erase event to keep the background transparent
    on windows.
    @param evt: wx.EVT_ERASE_BACKGROUND

    """
    pass

  def OnFocus(self, evt):
    """Set the visual focus state if need be"""
    if self._state['cur'] == GRADIENT_NORMAL:
        self.SetState(GRADIENT_HIGHLIGHT)

  def OnKeyUp(self, evt):
    """Execute a single button press action when the Return key is pressed
    and this control has the focus.
    @param evt: wx.EVT_KEY_UP

    """
    if evt.GetKeyCode() == wx.WXK_SPACE:
      self.SetState(GRADIENT_PRESSED)
      self.__PostEvent()
      wx.CallLater(100, self.SetState, GRADIENT_HIGHLIGHT)
    else:
      evt.Skip()

  def OnKillFocus(self, evt):
    """Set the visual state back to normal when focus is lost
    unless the control is currently in a pressed state.

    """
    # Note: this delay needs to be at least as much as the on in the KeyUp
    #       handler to prevent ghost highlighting from happening when
    #       quickly changing focus and activating buttons
    if self._state['cur'] != GRADIENT_PRESSED:
      self.SetState(GRADIENT_NORMAL)

  def OnLeftDown(self, evt):
    """Sets the pressed state and depending on the click position will
    show the popup menu if one has been set.

    """
    pos = evt.GetPositionTuple()
    self.SetState(GRADIENT_PRESSED)
    size = self.GetSizeTuple()
    if pos[0] >= size[0] - 16:
      if self._menu is not None:
        self.ShowMenu()

    self.SetFocus()

  def OnLeftUp(self, evt):
    """Post a button event if the control was previously in a
    pressed state.
    @param evt: wx.MouseEvent

    """
    if self._state['cur'] == GRADIENT_PRESSED:
      pos = evt.GetPositionTuple()
      size = self.GetSizeTuple()
      if self._disable_after_click > 0 :
        self.Enable(False)
      self.__PostEvent()
    self.SetState(GRADIENT_HIGHLIGHT)
    if self._disable_after_click > 0 :
      wx.CallLater(self._disable_after_click, lambda : self.Enable(True))

  def OnMenuClose(self, evt):
    """Refresh the control to a proper state after the menu has been
    dismissed.
    @param evt: wx.EVT_MENU_CLOSE

    """
    mpos = wx.GetMousePosition()
    if self.HitTest(self.ScreenToClient(mpos)) != wx.HT_WINDOW_OUTSIDE:
      self.SetState(GRADIENT_HIGHLIGHT)
    else:
      self.SetState(GRADIENT_NORMAL)
    evt.Skip()

  #---- End Event Handlers ----#

  def SetBitmap(self, bmp):
    """Set the bitmap displayed in the button
    @param bmp: wx.Bitmap

    """
    self._bmp['enable'] = bmp
    img = bmp.ConvertToImage()
    img = img.ConvertToGreyscale(.795, .073, .026) #(.634, .224, .143)
    self._bmp['disable'] = img.ConvertToBitmap()
    self.InvalidateBestSize()

  def SetBitmapDisabled(self, bmp):
    """Set the bitmap for the disabled state
    @param bmp: wx.Bitmap

    """
    self._bmp['disable'] = bmp

  # Aliases for SetBitmap* functions from BitmapButton
  SetBitmapFocus = SetBitmap
  SetBitmapHover = SetBitmap
  SetBitmapLabel = SetBitmap
  SetBitmapSelected = SetBitmap

  def SetFocus(self):
    """Set this control to have the focus"""
    if self._state['cur'] != GRADIENT_PRESSED:
      self.SetState(GRADIENT_HIGHLIGHT)
    wx.PyControl.SetFocus(self)

  def SetFont(self, font):
    """Adjust size of control when font changes"""
    wx.PyControl.SetFont(self, font)
    self.InvalidateBestSize()

  def SetLabel(self, label):
    """Set the label of the button
    @param label: lable string

    """
    wx.PyControl.SetLabel(self, label)
    self.InvalidateBestSize()

  def SetLabelColor(self, normal, hlight=wx.NullColour):
    """Set the color of the label. The optimal label color is usually
    automatically selected depending on the button color. In some
    cases the colors that are choosen may not be optimal.

    The normal state must be specified, if the other two params are left
    Null they will be automatically guessed based on the normal color. To
    prevent this automatic color choices from happening either specify
    a color or None for the other params.

    @param normal: Label color for normal state
    @keyword hlight: Color for when mouse is hovering over

    """
    self._color['default'] = False
    self.SetForegroundColour(normal)

    if hlight is not None:
      if hlight.IsOk():
        self._color['htxt'] = hlight
      else:
        self._color['htxt'] = BestLabelColour(normal)

    if wx.Platform == '__WXMSW__':
      self.GetParent().RefreshRect(self.GetRect(), False)
    else:
      self.Refresh()

  def SetMenu(self, menu):
    """Set the menu that can be shown when clicking on the
    drop arrow of the button.
    @param menu: wxMenu to use as a PopupMenu
    @note: Arrow is not drawn unless a menu is set

    """
    if self._menu is not None:
      self.Unbind(wx.EVT_MENU_CLOSE)

    self._menu = menu
    self.Bind(wx.EVT_MENU_CLOSE, self.OnMenuClose)
    self.InvalidateBestSize()

  def SetPressColor(self, color):
    """Set the color used for highlighting the pressed state
    @param color: wx.Color
    @note: also resets all text colours as necessary

    """
    self._color['default'] = False
    if color.Alpha() == 255:
      self._color['hlight'] = AdjustAlpha(color, 200)
    else:
      self._color['hlight'] = color
    #self._color['press'] = AdjustColour(color, -10, 160)
    self._color['htxt'] = BestLabelColour(self._color['hlight'])
    self.Refresh()

  def SetState(self, state):
    """Manually set the state of the button
    @param state: one of the MB_* values
    @note: the state may be altered by mouse actions

    """
    self._state['pre'] = self._state['cur']
    self._state['cur'] = state
    if wx.Platform == '__WXMSW__':
      self.GetParent().RefreshRect(self.GetRect(), False)
    else:
      self.Refresh()

  def SetWindowStyle(self, style):
    """Sets the window style bytes, the updates take place
    immediately no need to call refresh afterwards.
    @param style: bitmask of PB_STYLE_* values

    """
    self._style = style
    self.Refresh()

  def SetWindowVariant(self, variant):
    """Set the variant/font size of this control"""
    wx.PyControl.SetWindowVariant(self, variant)
    self.InvalidateBestSize()

  def ShouldInheritColours(self):
    """Overridden base class virtual. If the parent has non-default
    colours then we want this control to inherit them.

    """
    return True

  def ShowMenu(self):
    """Show the dropdown menu if one is associated with this control"""
    if self._menu is not None:
      size = self.GetSizeTuple()
      adj = wx.Platform == '__WXMAC__' and 3 or 0

      xpos = 1
      self.PopupMenu(self._menu, (xpos, size[1] + adj))

  def ToggleState(self):
    """Toggle button state"""
    if self._state['cur'] != GRADIENT_PRESSED:
      self.SetState(GRADIENT_PRESSED)
    else:
      self.SetState(GRADIENT_HIGHLIGHT)
