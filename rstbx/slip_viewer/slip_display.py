#!/usr/bin/env python
# -*- coding= utf-8 -*-

"""pySlip demonstration program."""
from __future__ import absolute_import, division, print_function
# Copyright (c) 2010, Ross Wilson (rzzzwilson@gmail.com). All rights reserved.
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import sys
import wx

from . import pyslip
from . import tile_generation
from six.moves import range
pyslip._Tiles = tile_generation._Tiles

######
# Various demo constants
######

# demo name/version
DemoName = 'pySlip %s - Demonstration' % pyslip.__version__
DemoVersion = '1.6'

# tiles info
TileDirectory = 'tiles'
MinTileLevel = 0

# initial view level and position
InitViewLevel = -2
#InitViewPosition = (103.770, 1.335)
InitViewPosition = (-100,100)

# levels on which various layers show
MRPointShowLevels = [3, 4]
MRImageShowLevels = [3, 4]
MRTextShowLevels = [3, 4]
MRPolyShowLevels = [3, 4]

# the number of decimal places in a lon/lat display
LonLatPrecision = 3

# startup size of the application
DefaultAppSize = (800, 665)

# how close click has to be before point is selected
# the value is distance squared (degrees^2)
PointSelectDelta = 0.025

# unselected point colour (rgb) and size
PointsColour = '#ff0000'
PointsSize = 3

# Selected point colour (rgb) and size
SelectColour = '#0000ff'
SelectSize = 5

# Polygon point colour (rgba) and size
PolygonColour = '#0000ff'
PolygonSize = 4

# Polygon2 point colour (rgba) and size
Polygon2Colour = '#000000'
Polygon2Size = 4

# image used for shipwrecks, glassy buttons, etc
ShipImg = 'graphics/shipwreck.png'

GlassyImg2 = 'graphics/glassy_button_2.png'
SelGlassyImg2 = 'graphics/selected_glassy_button_2.png'
GlassyImg3 = 'graphics/glassy_button_3.png'
SelGlassyImg3 = 'graphics/selected_glassy_button_3.png'
GlassyImg4 = 'graphics/glassy_button_4.png'
SelGlassyImg4 = 'graphics/selected_glassy_button_4.png'
GlassyImg5 = 'graphics/glassy_button_5.png'
SelGlassyImg5 = 'graphics/selected_glassy_button_5.png'
GlassyImg6 = 'graphics/glassy_button_6.png'
SelGlassyImg6 = 'graphics/selected_glassy_button_6.png'

# image used for shipwrecks
CompassRoseGraphic = 'graphics/compass_rose.png'

######
# Various GUI layout constants
######

# sizes of various spacers
HSpacerSize = (3,1)         # horizontal in application screen
VSpacerSize = (1,5)         # vertical in control pane

# border width when packing GUI elements
PackBorder = 1


################################################################################
# Override the wx.TextCtrl class to add read-only style and background colour
################################################################################

# background colour for the 'read-only' text field
ControlReadonlyColour = '#ffffcc'

class ROTextCtrl(wx.TextCtrl):
    """Override the wx.TextCtrl widget to get read-only text control which
    has a distinctive background colour."""

    def __init__(self, parent, value, tooltip='', *args, **kwargs):
        wx.TextCtrl.__init__(self, parent, wx.ID_ANY, value=value,
                             style=wx.TE_READONLY, *args, **kwargs)
        self.SetBackgroundColour(ControlReadonlyColour)
        self.SetToolTip(wx.ToolTip(tooltip))

################################################################################
# Override the wx.StaticBox class to show our style
################################################################################

class AppStaticBox(wx.StaticBox):

    def __init__(self, parent, label, *args, **kwargs):
        if label:
            label = '  ' + label + '  '
        wx.StaticBox.__init__(self, parent, wx.ID_ANY, label, *args, **kwargs)

################################################################################
# Class for a LayerControl widget.
################################################################################

myEVT_ONOFF = wx.NewEventType()
myEVT_SHOWONOFF = wx.NewEventType()
myEVT_SELECTONOFF = wx.NewEventType()

EVT_ONOFF = wx.PyEventBinder(myEVT_ONOFF, 1)
EVT_SHOWONOFF = wx.PyEventBinder(myEVT_SHOWONOFF, 1)
EVT_SELECTONOFF = wx.PyEventBinder(myEVT_SELECTONOFF, 1)

class LayerControlEvent(wx.PyCommandEvent):
    """Event sent when a LayerControl is changed."""

    def __init__(self, eventType, id):
        wx.PyCommandEvent.__init__(self, eventType, id)

class LayerControl(wx.Panel):

    def __init__(self, parent, title, selectable=False, editable=False,
                 **kwargs):
        """Initialise a LayerControl instance.

        parent      reference to parent object
        title       text to ahow in static box outline
        selectable  True if 'selectable' checkbox is to be displayed
        editable    True if layer can be edited
        **kwargs    keyword args for Panel
        """

        # create and initialise the base panel
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, **kwargs)
        self.SetBackgroundColour(wx.WHITE)

        self.selectable = selectable
        self.editable = editable

        box = AppStaticBox(self, title)
        sbs = wx.StaticBoxSizer(box, orient=wx.VERTICAL)
        gbs = wx.GridBagSizer()

        self.cbx_onoff = wx.CheckBox(self, wx.ID_ANY, label='Add layer')
        gbs.Add(self.cbx_onoff, (0,0), span=(1,4))

        self.cbx_show = wx.CheckBox(self, wx.ID_ANY, label='Show')
        gbs.Add(self.cbx_show, (1,1))
        self.cbx_show.Disable()

        if selectable:
            self.cbx_select = wx.CheckBox(self, wx.ID_ANY, label='Select')
            gbs.Add(self.cbx_select, (1,2))
            self.cbx_select.Disable()

        if editable:
            self.cbx_edit = wx.CheckBox(self, wx.ID_ANY, label='Edit')
            gbs.Add(self.cbx_edit, (1,3))
            self.cbx_edit.Disable()

        sbs.Add(gbs)
        self.SetSizer(sbs)
        sbs.Fit(self)

        # tie handlers to change events
        self.cbx_onoff.Bind(wx.EVT_CHECKBOX, self.onChangeOnOff)
        self.cbx_show.Bind(wx.EVT_CHECKBOX, self.onChangeShowOnOff)
        if selectable:
            self.cbx_select.Bind(wx.EVT_CHECKBOX, self.onChangeSelectOnOff)
#        if editable:
#            self.cbx_edit.Bind(wx.EVT_CHECKBOX, self.onChangeEditOnOff)

    def onChangeOnOff(self, event):
        """Main checkbox changed."""

        event = LayerControlEvent(myEVT_ONOFF, self.GetId())
        event.state = self.cbx_onoff.IsChecked()
        self.GetEventHandler().ProcessEvent(event)

        if self.cbx_onoff.IsChecked():
            self.cbx_show.Enable()
            self.cbx_show.SetValue(True)
            if self.selectable:
                self.cbx_select.Enable()
                self.cbx_select.SetValue(False)
            if self.editable:
                self.cbx_edit.Enable()
                self.cbx_edit.SetValue(False)
        else:
            self.cbx_show.Disable()
            if self.selectable:
                self.cbx_select.Disable()
            if self.editable:
                self.cbx_edit.Disable()

    def onChangeShowOnOff(self, event):
        """Show checkbox changed."""

        event = LayerControlEvent(myEVT_SHOWONOFF, self.GetId())
        event.state = self.cbx_show.IsChecked()
        self.GetEventHandler().ProcessEvent(event)

    def onChangeSelectOnOff(self, event):
        """Select checkbox changed."""

        event = LayerControlEvent(myEVT_SELECTONOFF, self.GetId())
        if self.selectable:
            event.state = self.cbx_select.IsChecked()
        else:
            event_state = False
        self.GetEventHandler().ProcessEvent(event)

################################################################################
# The main application frame
################################################################################

class AppFrame(wx.Frame):
    def __init__(self, tile_dir=TileDirectory):
        wx.Frame.__init__(self, None, size=DefaultAppSize,
                          title='%s %s' % (DemoName, DemoVersion))
        self.SetMinSize(DefaultAppSize)
        self.panel = wx.Panel(self, wx.ID_ANY)
        self.panel.SetBackgroundColour(wx.WHITE)
        self.panel.ClearBackground()

        self.tile_directory = tile_dir

        # build the GUI
        self.pyslip = None
        self.make_gui(self.panel)

        # do initialisation stuff - all the application stuff
        self.init()

        # finally, set up application window position
        self.Centre()

        # create select event dispatch directory
        self.demo_select_dispatch = {}

        # finally, bind events to handlers
        self.pyslip.Bind(pyslip.EVT_PYSLIP_SELECT, self.handle_select_event)
        self.pyslip.Bind(pyslip.EVT_PYSLIP_POSITION, self.handle_position_event)
        self.pyslip.Bind(pyslip.EVT_PYSLIP_LEVEL, self.handle_level_change)

#####
# Build the GUI
#####

    def set_pyslip(self, parent):
        self.pyslip = pyslip.PySlip(parent, tile_dir=None, # dummy file name value
                                    min_level=MinTileLevel)

    def make_gui(self, parent):
        """Create application GUI."""

        # start application layout
        all_display = wx.BoxSizer(wx.HORIZONTAL)
        parent.SetSizer(all_display)

        # put map view in left of horizontal box
        sl_box = self.make_gui_view(parent)
        all_display.Add(sl_box, proportion=1, border=1, flag=wx.EXPAND)

        # small spacer here - separate view and controls
        #all_display.AddSpacer(HSpacerSize)

        # add controls to right of spacer
        #controls = self.make_gui_controls(parent)
        #all_display.Add(controls, proportion=0, border=1)

        parent.SetSizerAndFit(all_display)

    def make_gui_view(self, parent):
        """Build the map view widget

        parent  reference to the widget parent

        Returns the static box sizer.
        """

        # create gui objects
        sb = AppStaticBox(parent, '')

        # lay out objects
        box = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        if (self.pyslip is None):
            self.set_pyslip(parent)
        box.Add(self.pyslip, proportion=1, border=1, flag=wx.EXPAND)

        return box

    def make_gui_controls(self, parent):
        """Build the 'controls' part of the GUI

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # all controls in vertical box sizer
        controls = wx.BoxSizer(wx.VERTICAL)

        # add the map level in use widget
        level = self.make_gui_level(parent)
        controls.Add(level, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # add the mouse position feedback stuff
        mouse = self.make_gui_mouse(parent)
        controls.Add(mouse, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for map-relative points layer
        point = self.make_gui_point(parent)
        controls.Add(point, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for view-relative points layer
        point_view = self.make_gui_point_view(parent)
        controls.Add(point_view, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for map-relative image layer
        image = self.make_gui_image(parent)
        controls.Add(image, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for map-relative image layer
        image_view = self.make_gui_image_view(parent)
        controls.Add(image_view, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for map-relative text layer
        text = self.make_gui_text(parent)
        controls.Add(text, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for view-relative text layer
        text_view = self.make_gui_text_view(parent)
        controls.Add(text_view, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for map-relative polygon layer
        poly = self.make_gui_poly(parent)
        controls.Add(poly, proportion=0, flag=wx.EXPAND|wx.ALL)

        # vertical spacer
        controls.AddSpacer(VSpacerSize)

        # controls for view-relative polygon layer
        poly_view = self.make_gui_poly_view(parent)
        controls.Add(poly_view, proportion=0, flag=wx.EXPAND|wx.ALL)

        return controls

    def make_gui_level(self, parent):
        """Build the control that shows the level.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create objects
        txt = wx.StaticText(parent, wx.ID_ANY, 'Level: ')
        self.map_level = wx.StaticText(parent, wx.ID_ANY, ' ')

        # lay out the controls
        sb = AppStaticBox(parent, 'Map level')
        box = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        box.Add(txt, border=PackBorder, flag=(wx.ALIGN_CENTER_VERTICAL
                                     |wx.ALIGN_RIGHT|wx.LEFT))
        box.Add(self.map_level, proportion=0, border=PackBorder,
                flag=wx.RIGHT|wx.TOP)

        return box

    def make_gui_mouse(self, parent):
        """Build the mouse part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create objects
        txt = wx.StaticText(parent, wx.ID_ANY, 'Lon/Lat: ')
        self.mouse_position = ROTextCtrl(parent, '', size=(150,-1),
                                         tooltip=('Shows the mouse '
                                                  'longitude and latitude '
                                                  'on the map'))

        # lay out the controls
        sb = AppStaticBox(parent, 'Mouse position')
        box = wx.StaticBoxSizer(sb, orient=wx.HORIZONTAL)
        box.Add(txt, border=PackBorder, flag=(wx.ALIGN_CENTER_VERTICAL
                                     |wx.ALIGN_RIGHT|wx.LEFT))
        box.Add(self.mouse_position, proportion=1, border=PackBorder,
                flag=wx.RIGHT|wx.TOP|wx.BOTTOM)

        return box

    def make_gui_point(self, parent):
        """Build the points part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        point_obj = LayerControl(parent, 'Points, map relative %s'
                                         % str(MRPointShowLevels),
                                 selectable=True)

        # tie to event handler(s)
        point_obj.Bind(EVT_ONOFF, self.pointOnOff)
        point_obj.Bind(EVT_SHOWONOFF, self.pointShowOnOff)
        point_obj.Bind(EVT_SELECTONOFF, self.pointSelectOnOff)

        return point_obj

    def make_gui_point_view(self, parent):
        """Build the view-relative points part of the GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        point_obj = LayerControl(parent, 'Points, view relative',
                                 selectable=True)

        # tie to event handler(s)
        point_obj.Bind(EVT_ONOFF, self.pointViewOnOff)
        point_obj.Bind(EVT_SHOWONOFF, self.pointViewShowOnOff)
        point_obj.Bind(EVT_SELECTONOFF, self.pointViewSelectOnOff)

        return point_obj

    def make_gui_image(self, parent):
        """Build the image part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        image_obj = LayerControl(parent, 'Images, map relative %s'
                                         % str(MRImageShowLevels),
                                 selectable=True)

        # tie to event handler(s)
        image_obj.Bind(EVT_ONOFF, self.imageOnOff)
        image_obj.Bind(EVT_SHOWONOFF, self.imageShowOnOff)
        image_obj.Bind(EVT_SELECTONOFF, self.imageSelectOnOff)

        return image_obj

    def make_gui_image_view(self, parent):
        """Build the view-relative image part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        image_obj = LayerControl(parent, 'Images, view relative')

        # tie to event handler(s)
        image_obj.Bind(EVT_ONOFF, self.imageViewOnOff)
        image_obj.Bind(EVT_SHOWONOFF, self.imageViewShowOnOff)
        image_obj.Bind(EVT_SELECTONOFF, self.imageViewSelectOnOff)

        return image_obj

    def make_gui_text(self, parent):
        """Build the map-relative text part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        text_obj = LayerControl(parent, 'Text, map relative %s'
                                        % str(MRTextShowLevels),
                                selectable=True, editable=True)

        # tie to event handler(s)
        text_obj.Bind(EVT_ONOFF, self.textOnOff)
        text_obj.Bind(EVT_SHOWONOFF, self.textShowOnOff)
        text_obj.Bind(EVT_SELECTONOFF, self.textSelectOnOff)

        return text_obj

    def make_gui_text_view(self, parent):
        """Build the view-relative text part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        text_view_obj = LayerControl(parent, 'Text, view relative',
                                     selectable=True)

        # tie to event handler(s)
        text_view_obj.Bind(EVT_ONOFF, self.textViewOnOff)
        text_view_obj.Bind(EVT_SHOWONOFF, self.textViewShowOnOff)
        text_view_obj.Bind(EVT_SELECTONOFF, self.textViewSelectOnOff)

        return text_view_obj

    def make_gui_poly(self, parent):
        """Build the map-relative polygon part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        poly_obj = LayerControl(parent,
                                'Polygon, map relative %s'
                                     % str(MRPolyShowLevels),
                                selectable=True)

        # tie to event handler(s)
        poly_obj.Bind(EVT_ONOFF, self.polyOnOff)
        poly_obj.Bind(EVT_SHOWONOFF, self.polyShowOnOff)
        poly_obj.Bind(EVT_SELECTONOFF, self.polySelectOnOff)

        return poly_obj

    def make_gui_poly_view(self, parent):
        """Build the view-relative polygon part of the controls part of GUI.

        parent  reference to parent

        Returns reference to containing sizer object.
        """

        # create widgets
        poly_view_obj = LayerControl(parent, 'Polygon, view relative',
                                     selectable=True)

        # tie to event handler(s)
        poly_view_obj.Bind(EVT_ONOFF, self.polyViewOnOff)
        poly_view_obj.Bind(EVT_SHOWONOFF, self.polyViewShowOnOff)
        poly_view_obj.Bind(EVT_SELECTONOFF, self.polyViewSelectOnOff)

        return poly_view_obj

    ######
    # pySlip demo control event handlers
    ######

##### map-relative point layer

    def pointOnOff(self, event):
        """Handle OnOff event for point layer control."""

        if event.state:
            self.point_layer = \
                self.pyslip.AddPointLayer(PointData, map_rel=True,
                                          color=PointDataColour, radius=3,
                                          offset_x=0, offset_y=0, visible=True,
                                          show_levels=MRPointShowLevels,
                                          name='<pt_layer>')
        else:
            self.pyslip.DeleteLayer(self.point_layer)
            self.point_layer = None
            if self.sel_point_layer:
                self.pyslip.DeleteLayer(self.sel_point_layer)
                self.sel_point_layer = None
                self.sel_point = None

    def pointShowOnOff(self, event):
        """Handle ShowOnOff event for point layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.point_layer)
            if self.sel_point_layer:
                self.pyslip.ShowLayer(self.sel_point_layer)
        else:
            self.pyslip.HideLayer(self.point_layer)
            if self.sel_point_layer:
                self.pyslip.HideLayer(self.sel_point_layer)

    def pointSelectOnOff(self, event):
        """Handle SelectOnOff event for point layer control."""

        layer = self.point_layer
        if event.state:
            self.add_select_handler(layer, self.pointSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def pointSelect(self, event):
        """Handle point select exception from pyslip.

        event  the event that contains these attributes:
                   layer_id  ID of the layer the select is for
                   sel_type  type of select event
                   point     selected point(s) geo coords+data
                                 ((x,y), data)
                             (if None then no point(s) selected)

        The point select is designed to be click for on,
        then click again for off.
        """

        if event.evtype == pyslip.EventPointSelect:
            if event.point:
                (point, data) = event.point
            if event.point is None or point == self.sel_point:
                # select again, turn point off
                self.sel_point = None
                self.pyslip.DeleteLayer(self.sel_point_layer)
                self.sel_point_layer = None
            elif point:
                if self.sel_point_layer:
                    self.pyslip.DeleteLayer(self.sel_point_layer)
                self.sel_point = point
                self.sel_point_layer = \
                    self.pyslip.AddPointLayer((point,), map_rel=True,
                                              color='#0000ff',
                                              radius=5, visible=True,
                                              show_levels=MRPointShowLevels,
                                              name='<sel_pt_layer>')
        if event.evtype == pyslip.EventBoxSelect: # left box select
            # remove any previous selection
            if self.sel_point_layer:
                self.pyslip.DeleteLayer(self.sel_point_layer)
                self.sel_point_layer = None

            if event.point:
                pts = [pt for (pt,d) in event.point]
                self.sel_point_layer = \
                    self.pyslip.AddPointLayer(pts, map_rel=True,
                                              color='#00ffff',
                                              radius=5, visible=True,
                                              show_levels=[3,4],
                                              name='<boxsel_pt_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_point_layer,
                                                 self.point_layer)

        return True

##### view-relative point layer

    def pointViewOnOff(self, event):
        """Handle OnOff event for point view layer control."""

        if event.state:
            self.point_view_layer = \
                self.pyslip.AddPointLayer(PointViewData, map_rel=False,
                                          placement='se',
                                          color=PointViewDataColour, radius=1,
                                          visible=True,
                                          name='<point_view_layer>')
        else:
            self.pyslip.DeleteLayer(self.point_view_layer)
            self.point_view_layer = None
            if self.sel_point_view_layer:
                self.pyslip.DeleteLayer(self.sel_point_view_layer)
                self.sel_point_view_layer = None
                self.sel_point_view = None

    def pointViewShowOnOff(self, event):
        """Handle ShowOnOff event for point view layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.point_view_layer)
            if self.sel_point_view_layer:
                self.pyslip.ShowLayer(self.sel_point_view_layer)
        else:
            self.pyslip.HideLayer(self.point_view_layer)
            if self.sel_point_view_layer:
                self.pyslip.HideLayer(self.sel_point_view_layer)

    def pointViewSelectOnOff(self, event):
        """Handle SelectOnOff event for point view layer control."""

        layer = self.point_view_layer
        if event.state:
            self.add_select_handler(layer, self.pointViewSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def pointViewSelect(self, event):
        """Handle view-relative point select event from pyslip.

        event  the event that contains these attributes:
                   layer_id  ID of the layer the select is for
                   sel_type  type of select event
                   point     selected point(s) geo coordinates
                             (if None then no point(s) was selected)

        The point select is designed to be click for on,
        then click again for off.
        """


        if event.evtype == pyslip.EventPointSelect:
            if event.point:
                (point, data) = event.point
            if event.point is None or point == self.sel_point_view:
                # select again, turn point off
                self.sel_point_view = None
                self.pyslip.DeleteLayer(self.sel_point_view_layer)
                self.sel_point_view_layer = None
            elif event.point:
                if self.sel_point_view_layer:
                    self.pyslip.DeleteLayer(self.sel_point_view_layer)
                self.sel_point_view = point
                self.sel_point_view_layer = \
                    self.pyslip.AddPointLayer((point,), map_rel=False,
                                              color='#0000ff',
                                              radius=3, visible=True,
                                              name='<sel_pt_view_layer>')
        elif event.evtype == pyslip.EventBoxSelect:
            # remove any previous selection
            if self.sel_point_view_layer:
                self.pyslip.DeleteLayer(self.sel_point_view_layer)
                self.sel_point_view_layer = None

            if event.point:
                pts = [pt for (pt,d) in event.point]
                self.sel_point_view_layer = \
                    self.pyslip.AddPointLayer(pts, map_rel=False,
                                              color='#00ffff',
                                              radius=3, visible=True,
                                              name='<boxsel_pt_view_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_point_view_layer,
                                                 self.point_view_layer)

        return True

##### map-relative image layer

    def imageOnOff(self, event):
        """Handle OnOff event for map-relative image layer control."""

        if event.state:
            self.image_layer = \
                self.pyslip.AddImageLayer(ImageData, map_rel=True,
                                          visible=True,
                                          show_levels=MRImageShowLevels,
                                          name='<image_layer>')
        else:
            self.pyslip.DeleteLayer(self.image_layer)
            self.image_layer = None
            if self.sel_image_layer:
                self.pyslip.DeleteLayer(self.sel_image_layer)
                self.sel_image_layer = None
                self.sel_image = None

    def imageShowOnOff(self, event):
        """Handle ShowOnOff event for image layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.image_layer)
            if self.sel_image_layer:
                self.pyslip.ShowLayer(self.sel_image_layer)
        else:
            self.pyslip.HideLayer(self.image_layer)
            if self.sel_image_layer:
                self.pyslip.HideLayer(self.sel_image_layer)

    def imageSelectOnOff(self, event):
        """Handle SelectOnOff event for image layer control."""

        layer = self.image_layer
        if event.state:
            self.add_select_handler(layer, self.imageSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def imageSelect(self, event):
        """Select event from pyslip."""

        point = event.point

        if event.evtype == pyslip.EventPointSelect:
            if point == self.sel_image:
                # select again, turn point off
                self.sel_image = None
                self.pyslip.DeleteLayer(self.sel_image_layer)
                self.sel_image_layer = None
            elif point:
                if self.sel_image_layer:
                    self.pyslip.DeleteLayer(self.sel_image_layer)
                self.sel_image = point
                self.sel_image_layer = \
                    self.pyslip.AddPointLayer((point,), map_rel=True,
                                              color='#0000ff',
                                              radius=5, visible=True,
                                              show_levels=[3,4],
                                              name='<sel_pt_layer>')
        elif event.evtype == pyslip.EventBoxSelect:
            # remove any previous selection
            if self.sel_image_layer:
                self.pyslip.DeleteLayer(self.sel_image_layer)
                self.sel_image_layer = None

            if point:
                self.sel_image_layer = \
                    self.pyslip.AddPointLayer(point, map_rel=True,
                                              color='#00ffff',
                                              radius=5, visible=True,
                                              show_levels=[3,4],
                                              name='<boxsel_pt_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_image_layer,
                                                 self.image_layer)

        return True

    def imageBSelect(self, id, points=None):
        """Select event from pyslip."""

        # remove any previous selection
        if self.sel_image_layer:
            self.pyslip.DeleteLayer(self.sel_image_layer)
            self.sel_image_layer = None

        if points:
            self.sel_image_layer = \
                self.pyslip.AddPointLayer(points, map_rel=True,
                                          color='#e0e0e0',
                                          radius=13, visible=True,
                                          show_levels=[3,4],
                                          name='<boxsel_img_layer>')
            self.pyslip.PlaceLayerBelowLayer(self.sel_image_layer,
                                             self.image_layer)

        return True

##### view-relative image layer

    def imageViewOnOff(self, event):
        """Handle OnOff event for view-relative image layer control."""

        if event.state:
            self.image_view_layer = \
                self.pyslip.AddImageLayer(ImageViewData, map_rel=False,
                                          visible=True,
                                          name='<image_view_layer>')
        else:
            self.pyslip.DeleteLayer(self.image_view_layer)
            self.image_view_layer = None
            if self.sel_image_view_layer:
                self.pyslip.DeleteLayer(self.sel_image_view_layer)
                self.sel_image_view_layer = None
                self.sel_image_view_point = None

    def imageViewShowOnOff(self, event):
        """Handle ShowOnOff event for image layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.image_view_layer)
            if self.sel_image_view_layer:
                self.pyslip.ShowLayer(self.sel_image_layer)
        else:
            self.pyslip.HideLayer(self.image_view_layer)
            if self.sel_image_view_layer:
                self.pyslip.HideLayer(self.sel_image_layer)

    def imageViewSelectOnOff(self, event):
        """Handle SelectOnOff event for image layer control."""

        layer = self.image_view_layer
        if event.state:
            self.add_select_handler(layer, self.imageViewSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def imageViewSelect(self, id, posn=None):
        """View-relative image select event from pyslip."""

        if posn:
            for p in ImageViewData:
                pp = (p[0], p[1])
                if pp == posn:
                    if pp == self.sel_image_view:
                        # select again, turn point off
                        self.sel_image_view = None
                        self.pyslip.DeleteLayer(self.sel_image_view_layer)
                        self.sel_image_view_layer = None
                    else:
                        if self.sel_image_view_layer:
                            self.pyslip.DeleteLayer(self.sel_image_view_layer)
                        self.sel_image_view = pp
                        self.sel_image_view_layer = \
                            self.pyslip.AddPointLayer((pp,), map_rel=False,
                                                      color='#00ffff',
                                                      radius=5, visible=True,
                                                      name='<sel_image_view>')
        return True

##### map-relative text layer

    def textOnOff(self, event):
        """Handle OnOff event for map-relative text layer control."""

        if event.state:
            self.text_layer = \
                self.pyslip.AddTextLayer(TextData, map_rel=True,
                                         name='<text_layer>', visible=True,
                                         show_levels=MRTextShowLevels,
                                         placement='ne')
        else:
            self.pyslip.DeleteLayer(self.text_layer)
            self.text_layer = None
            if self.sel_text_layer:
                self.pyslip.DeleteLayer(self.sel_text_layer)
                self.sel_text_layer = None
                self.sel_text_point = None

    def textShowOnOff(self, event):
        """Handle ShowOnOff event for text layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.text_layer)
            if self.sel_text_layer:
                self.pyslip.ShowLayer(self.sel_text_layer)
        else:
            self.pyslip.HideLayer(self.text_layer)
            if self.sel_text_layer:
                self.pyslip.HideLayer(self.sel_text_layer)

    def textSelectOnOff(self, event):
        """Handle SelectOnOff event for text layer control."""

        layer = self.text_layer
        if event.state:
            self.add_select_handler(layer, self.textSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)


    def textSelect(self, event):
        """Map-relative text select event from pyslip."""

        if event.evtype == pyslip.EventPointSelect:
            if event.point:
                (point, data) = event.point
            if event.point is None or point == self.sel_text:
                # select again, turn point off
                self.sel_text = None
                self.pyslip.DeleteLayer(self.sel_text_layer)
                self.sel_text_layer = None
            elif point:
                if self.sel_text_layer:
                    self.pyslip.DeleteLayer(self.sel_text_layer)
                self.sel_text = point
                self.sel_text_layer = \
                    self.pyslip.AddPointLayer((point,), map_rel=True,
                                              color='#0000ff',
                                              radius=5, visible=True,
                                              show_levels=MRTextShowLevels,
                                              name='<sel_text_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_text_layer,
                                                 self.text_layer)
        if event.evtype == pyslip.EventBoxSelect: # left box select
            # remove any previous selection
            if self.sel_text_layer:
                self.pyslip.DeleteLayer(self.sel_text_layer)
                self.sel_text_layer = None

            if event.point:
                pts = [pt for (pt,d) in event.point]
                self.sel_text_layer = \
                    self.pyslip.AddPointLayer(pts, map_rel=True,
                                              color='#00ffff',
                                              radius=5, visible=True,
                                              show_levels=[3,4],
                                              name='<boxsel_text_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_text_layer,
                                                 self.text_layer)

        return True

##### view-relative text layer

    def textViewOnOff(self, event):
        """Handle OnOff event for map-relative text layer control."""

        if event.state:
            self.text_view_layer = \
                self.pyslip.AddTextLayer(TextViewData, map_rel=False,
                                         name='<text_view_layer>',
                                         placement='cn', visible=True,
                                         fontsize=24, textcolor='#0000ff')
        else:
            self.pyslip.DeleteLayer(self.text_view_layer)
            self.text_view_layer = None
            if self.sel_text_view_layer:
                self.pyslip.DeleteLayer(self.sel_text_view_layer)
                self.sel_text_view_layer = None
                self.sel_text_view_point = None

    def textViewShowOnOff(self, event):
        """Handle ShowOnOff event for text layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.text_view_layer)
            if self.sel_text_view_layer:
                self.pyslip.ShowLayer(self.sel_text_view_layer)
        else:
            self.pyslip.HideLayer(self.text_view_layer)
            if self.sel_text_view_layer:
                self.pyslip.HideLayer(self.sel_text_view_layer)

    def textViewSelectOnOff(self, event):
        """Handle SelectOnOff event for text layer control."""

        layer = self.text_view_layer
        if event.state:
            self.add_select_handler(layer, self.textViewSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def textViewSelect(self, event):
        """Map-relative text select event from pyslip."""

#        if posn:
#            for p in TextData:
#                pp = (p[0], p[1])
#                if pp == posn:
#                    if pp == self.sel_view_text:
#                        # select again, turn point off
#                        self.sel_view_text = None
#                        self.pyslip.DeleteLayer(self.sel_text_view_layer)
#                        self.sel_text_view_layer = None
#                    else:
#                        if self.sel_text_view_layer:
#                            self.pyslip.DeleteLayer(self.sel_text_view_layer)
#                        self.sel_view_text = pp
#                        self.sel_text_view_layer = \
#                            self.pyslip.AddPointLayer((pp,), map_rel=True,
#                                                      color='#80ffff',
#                                                      radius=5, visible=True,
#                                                      name='<sel_text>')
#        return True

        if event.evtype == pyslip.EventPointSelect:
            if event.point:
                (point, data) = event.point
            if event.point is None or point == self.sel_text_view:
                # select again, turn point off
                self.sel_text_view = None
                self.pyslip.DeleteLayer(self.sel_text_view_layer)
                self.sel_text_view_layer = None
            elif point:
                if self.sel_text_view_layer:
                    self.pyslip.DeleteLayer(self.sel_text_view_layer)
                self.sel_text_view = point
                self.sel_text_view_layer = \
                    self.pyslip.AddPointLayer((point,), map_rel=True,
                                              color='#80ffff',
                                              radius=5, visible=True,
                                              show_levels=MRTextShowLevels,
                                              name='<sel_text_view_layer>')
        elif event.evtype == pyslip.EventRightPointSelect:    # right point select
            pass
        if event.evtype == pyslip.EventBoxSelect: # left box select
            # remove any previous selection
            if self.sel_text_view_layer:
                self.pyslip.DeleteLayer(self.sel_text_view_layer)
                self.sel_text_view_layer = None

            if event.point:
                pts = [pt for (pt,d) in event.point]
                self.sel_text_view_layer = \
                    self.pyslip.AddPointLayer(pts, map_rel=True,
                                              color='#00ffff',
                                              radius=5, visible=True,
                                              show_levels=[3,4],
                                              name='<boxsel_text_view_layer>')
                self.pyslip.PlaceLayerBelowLayer(self.sel_text_view_layer,
                                                 self.text_view_layer)
        else:   # right box select
            pass

        return True

##### map-relative polygon layer

    def polyOnOff(self, event):
        """Handle OnOff event for map-relative polygon layer control."""

        if event.state:
            self.poly_layer = \
                self.pyslip.AddPolygonLayer(PolyData, map_rel=True,
                                            visible=True,
                                            show_levels=MRPolyShowLevels,
                                            name='<poly_layer>')
        else:
            self.pyslip.DeleteLayer(self.poly_layer)
            self.poly_layer = None
            if self.sel_poly_layer:
                self.pyslip.DeleteLayer(self.sel_poly_layer)
                self.sel_poly_layer = None
                self.sel_poly_point = None

    def polyShowOnOff(self, event):
        """Handle ShowOnOff event for polygon layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.poly_layer)
            if self.sel_poly_layer:
                self.pyslip.ShowLayer(self.sel_poly_layer)
        else:
            self.pyslip.HideLayer(self.poly_layer)
            if self.sel_poly_layer:
                self.pyslip.HideLayer(self.sel_poly_layer)

    def polySelectOnOff(self, event):
        """Handle SelectOnOff event for polygon layer control."""

        layer = self.poly_layer
        if event.state:
            self.add_select_handler(layer, self.polySelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def polySelect(self, event):
        """Map-relative polygon select event from pyslip.

        event  the pyslip event which has attributes:
                   evtype    the event type
                   layer_id  ID of the layer selected
                   point     iterable of poly points, can be None
                   mposn     map-relative position of mouse-click
                   vposn     view-relative position of mouse-click
        """

        poly = event.point

        if event.evtype == pyslip.EventPointSelect:
            if poly:
                if poly == self.sel_poly:
                    # polygon selected again, turn selection off
                    self.sel_poly = None
                    self.pyslip.DeleteLayer(self.sel_poly_layer)
                    self.sel_poly_layer = None
                else:
                    # new selection
                    if self.sel_poly_layer:
                        # deselect previously selected poly
                        self.pyslip.DeleteLayer(self.sel_poly_layer)
                    self.sel_poly = poly
                    self.sel_poly_layer = \
                        self.pyslip.AddPointLayer(poly, map_rel=True,
                                                  color='#ff00ff',
                                                  radius=9, visible=True,
                                                  show_levels=[3,4],
                                                  name='<sel_poly>')
        else:   # box select, not yet implemented
            pass

        return True

##### view-relative polygon layer

    def polyViewOnOff(self, event):
        """Handle OnOff event for map-relative polygon layer control."""

        if event.state:
            self.poly_view_layer = \
                self.pyslip.AddPolygonLayer(PolyViewData, map_rel=False,
                                            name='<poly_view_layer>',
                                            placement='cn', visible=True,
                                            fontsize=24, color='#0000ff')
        else:
            self.pyslip.DeleteLayer(self.poly_view_layer)
            self.poly_view_layer = None
            if self.sel_poly_view_layer:
                self.pyslip.DeleteLayer(self.sel_poly_view_layer)
                self.sel_poly_view_layer = None
                self.sel_poly_view_point = None

    def polyViewShowOnOff(self, event):
        """Handle ShowOnOff event for polygon layer control."""

        if event.state:
            self.pyslip.ShowLayer(self.poly_view_layer)
            if self.sel_poly_view_layer:
                self.pyslip.ShowLayer(self.sel_poly_view_layer)
        else:
            self.pyslip.HideLayer(self.poly_view_layer)
            if self.sel_poly_view_layer:
                self.pyslip.HideLayer(self.sel_poly_view_layer)

    def polyViewSelectOnOff(self, event):
        """Handle SelectOnOff event for polygon layer control."""

        layer = self.poly_view_layer
        if event.state:
            self.add_select_handler(layer, self.polyViewSelect)
            self.pyslip.SetLayerSelectable(layer, True)
        else:
            self.del_select_handler(layer)
            self.pyslip.SetLayerSelectable(layer, False)

    def polyViewSelect(self, id, posn=None):
        """View-relative polygon select event from pyslip."""

        if posn:
            for p in PolyData:
                pp = (p[0], p[1])
                if pp == posn:
                    if pp == self.sel_view_poly:
                        # select again, turn polygon off
                        self.view_sel_poly = None
                        self.pyslip.DeleteLayer(self.sel_poly_view_layer)
                        self.sel_poly_view_layer = None
                    else:
                        if self.sel_poly_view_layer:
                            self.pyslip.DeleteLayer(self.sel_poly_view_layer)
                        self.view_sel_poly = pp
                        self.sel_poly_view_layer = \
                            self.pyslip.AddPointLayer((pp,), map_rel=True,
                                                      color='#80ffff',
                                                      radius=5, visible=True,
                                                      name='<sel_polyv>')
        return True

    ######
    # Finish initialization of data, etc
    ######
    def init_image_viewer(self):
        self.pyslip.GotoLevelAndPosition(InitViewLevel, InitViewPosition)


    def init(self):
        global PointData, PointDataColour
        global PointViewData, PointViewDataColour
        global ImageData
        global ImageViewData
        global TextData # , TextDataColour
        global TextViewData
        global PolyData
        global PolyViewData

        # create PointData
        PointData = []
        count = 0
        for lon in range(-70, 290+1, 5):
            for lat in range(-65, 65+1, 5):
                PointData.append((lon, lat, {'data': count}))
                count += 1
        PointDataColour = '#ff000080'   # semi-transparent

        # create PointViewData - a point-rendition of 'PYSLIP'
        PointViewData = [(-66,-14),(-66,-13),(-66,-12),(-66,-11),(-66,-10),
                         (-66,-9),(-66,-8),(-66,-7),(-66,-6),(-66,-5),(-66,-4),
                         (-66,-3),(-65,-7),(-64,-7),(-63,-7),(-62,-7),(-61,-8),
                         (-60,-9),(-60,-10),(-60,-11),(-60,-12),(-61,-13),
                         (-62,-14),(-63,-14),(-64,-14),(65,-14),            # P
                         (-59,-14),(-58,-13),(-57,-12),(-56,-11),(-55,-10),
                         (-53,-10),(-52,-11),(-51,-12),(-50,-13),(-49,-14),
                         (-54,-9),(-54,-8),(-54,-7),(-54,-6),(-54,-5),
                         (-54,-4),(-54,-3),                                 # Y
                         (-41,-13),(-42,-14),(-43,-14),(-44,-14),(-45,-14),
                         (-46,-14),(-47,-13),(-48,-12),(-48,-11),(-47,-10),
                         (-46,-9),(-45,-9),(-44,-9),(-43,-9),(-42,-8),
                         (-41,-7),(-41,-6),(-41,-5),(-42,-4),(-43,-3),
                         (-44,-3),(-45,-3),(-46,-3),(-47,-3),(-48,-4),       # S
                         (-39,-14),(-39,-13),(-39,-12),(-39,-11),(-39,-10),
                         (-39,-9),(-39,-8),(-39,-7),(-39,-6),(-39,-5),
                         (-39,-4),(-39,-3),(-38,-3),(-37,-3),(-36,-3),
                         (-35,-3),(-34,-3),(-33,-3),(-32,-3),                # L
                         (-29,-14),(-29,-13),(-29,-12),
                         (-29,-11),(-29,-10),(-29,-9),(-29,-8),(-29,-7),
                         (-29,-6),(-29,-5),(-29,-4),(-29,-3),                # I
                         (-26,-14),(-26,-13),(-26,-12),(-26,-11),(-26,-10),
                         (-26,-9),(-26,-8),(-26,-7),(-26,-6),(-26,-5),(-26,-4),
                         (-26,-3),(-25,-7),(-24,-7),(-23,-7),(-22,-7),(-21,-8),
                         (-20,-9),(-20,-10),(-20,-11),(-20,-12),(-21,-13),
                         (-22,-14),(-23,-14),(-24,-14),(25,-14)]             # P
        PointViewDataColour = '#00ff0020'       # very transparent

        # create image data
        ImageData = [# Agnes Napier - 1855
                     (160.0, -30.0, ShipImg, {'placement': 'cc'}),
                     # Venus - 1826
                     (145.0, -11.0, ShipImg, {'placement': 'ne'}),
                     # Wolverine - 1879
                     (156.0, -23.0, ShipImg, {'placement': 'nw'}),
                     # Thomas Day - 1884
                     (150.0, -15.0, ShipImg, {'placement': 'sw'}),
                     # Sybil - 1902
                     (165.0, -19.0, ShipImg, {'placement': 'se'}),
                     # Prince of Denmark - 1863
                     (158.55, -19.98, ShipImg),
                     # Moltke - 1911
                     (146.867525, -19.152185, ShipImg)
                    ]
        ImageData2 = []
        ImageData3 = []
        ImageData4 = []
        ImageData5 = []
        ImageData6 = []
        self.map_level_2_img = {0: ImageData2,
                                1: ImageData3,
                                2: ImageData4,
                                3: ImageData5,
                                4: ImageData6}
        self.map_level_2_selimg = {0: SelGlassyImg2,
                                   1: SelGlassyImg3,
                                   2: SelGlassyImg4,
                                   3: SelGlassyImg5,
                                   4: SelGlassyImg6}
        self.current_layer_img_layer = None
        for x in range(80):
            for y in range(40):
                ImageData.append((-30+x*2, y*2-30, GlassyImg4))

        ImageViewData = [(0, 0, CompassRoseGraphic, {'placement': 'ne'})]

        text_placement = {'placement': 'se'}
        transparent_placement = {'placement': 'se', 'colour': '#00000040'}
        capital = {'placement': 'se', 'fontsize': 14, 'color': 'red',
                   'textcolour': 'red'}
        TextData = [(151.20, -33.85, 'Sydney', text_placement),
                    (144.95, -37.84, 'Melbourne', {'placement': 'ce'}),
                    (153.08, -27.48, 'Brisbane', text_placement),
                    (115.86, -31.96, 'Perth', transparent_placement),
                    (138.30, -35.52, 'Adelaide', text_placement),
                    (130.98, -12.61, 'Darwin', text_placement),
                    (147.31, -42.96, 'Hobart', text_placement),
                    (174.75, -36.80, 'Auckland', text_placement),
                    (174.75, -41.29, 'Wellington', capital),
                    (172.61, -43.51, 'Christchurch', text_placement),
                    (168.74, -45.01, 'Queenstown', text_placement),
                    (147.30, -09.41, 'Port Moresby', capital),
                    (106.822922, -6.185451, 'Jakarta', capital),
                    (110.364444, -7.801389, 'Yogyakarta', text_placement),
                    (120.966667, 14.563333, 'Manila', capital),
                    (271.74, +40.11, 'Champaign', text_placement),
                    (160.0, -30.0, 'Agnes Napier - 1855',
                        {'placement': 'cw', 'offset_x': 20, 'color': 'green'}),
                    (145.0, -11.0, 'Venus - 1826',
                        {'placement': 'sw', 'color': 'green'}),
                    (156.0, -23.0, 'Wolverine - 1879',
                        {'placement': 'ce', 'color': 'green'}),
                    (150.0, -15.0, 'Thomas Day - 1884',
                        {'color': 'green'}),
                    (165.0, -19.0, 'Sybil - 1902',
                        {'placement': 'cw', 'color': 'green'}),
                    (158.55, -19.98, 'Prince of Denmark - 1863',
                        {'placement': 'nw', 'offset_x': 20, 'color': 'green'}),
                    (146.867525, -19.152182, 'Moltke - 1911',
                        {'placement': 'ce', 'offset_x': 20, 'color': 'green'})
                   ]
        if sys.platform != 'win32':
            TextData.extend([
                    (106.36, +10.36, 'Mỹ Tho', {'placement': 'ne'}),
                    (105.85, +21.033333, 'Hà Nội', capital),
                    (106.681944, 10.769444, 'Thành phố Hồ Chí Minh', {'placement': 'sw'}),
                    (132.47, +34.44, '広島市 (Hiroshima City)', text_placement),
                    (114.158889, +22.278333, '香港 (Hong Kong)', {'placement': 'nw'}),
                    ( 96.16, +16.80, 'ရန်ကုန် (Yangon)', capital),
                    (104.93, +11.54, ' ភ្នំពេញ (Phnom Penh)',
                        {'placement': 'ce', 'fontsize': 12, 'color': 'red'}),
                    (100.49, +13.75, 'กรุงเทพมหานคร (Bangkok)', capital),
                    ( 77.56, +34.09, 'གླེ་(Leh)', text_placement),
                    (84.991275, 24.695102, 'बोधगया (Bodh Gaya)', text_placement)])
#        TextDataColour = '#ffffff40'

        TextViewData = [(0, 7, '%s %s' % (DemoName, DemoVersion))]

        PolyData = [(((150,10),(160,20),(170,10),(165,0),(155,0)),
                      {'width': 3, 'color': 'blue', 'closed': True}),
                    (((165,-35),(175,-35),(175,-45),(165,-45)),
                      {'width': 10, 'color': '#00ff00c0', 'filled': True,
                       'fillcolor': '#ffff0040'}),
                    (((190,-30),(220,-50),(220,-30),(190,-50)),
                      {'width': 3, 'color': 'green', 'filled': True,
                       'fillcolor': 'yellow'}),
                    (((190,+50),(220,+65),(220,+50),(190,+65)),
                      {'width': 10, 'color': '#00000040'})
                   ]

        PolyViewData = [(((0,0),(230,0),(230,40),(-230,40),(-230,0)),
                        {'width': 3, 'color': '#00ff00ff', 'closed': True,
                         'placement': 'cn', 'offset_y': 1})]

        # set initial view position
        self.pyslip.GotoLevelAndPosition(InitViewLevel, InitViewPosition)
        self.map_level.SetLabel('%d' % InitViewLevel)

        # define layer ID variables & sub-checkbox state variables
        self.point_layer = None
        self.sel_point_layer = None
        self.sel_point = None

        self.point_view_layer = None
        self.sel_point_view_layer = None
        self.sel_point_view = None

        self.image_layer = None
        self.sel_image_layer = None
        self.sel_image = None

        self.image_view_layer = None
        self.sel_image_view_layer = None
        self.sel_image_view = None

        self.text_layer = None
        self.sel_text_layer = None
        self.sel_text = None

        self.text_view_layer = None
        self.sel_text_view_layer = None
        self.sel_view_text = None

        self.poly_layer = None
        self.sel_poly_layer = None
        self.sel_poly = None

        self.poly_view_layer = None
        self.sel_poly_view_layer = None
        self.sel_poly = None

        # force pyslip initialisation
        self.pyslip.OnSize()

    ######
    # Exception handlers
    ######

    def handle_select_event(self, event):
        """Handle a pySlip point/box SELECT event."""

        layer_id = event.layer_id
        #point = event.point

        self.demo_select_dispatch.get(layer_id, self.null_handler)(event)

    def null_handler(self, event):
        """Routine to handle unexpected events."""

        print('ERROR: null_handler!?')

    def handle_position_event(self, event):
        """Handle a pySlip POSITION event."""

        posn_str = ''
        if event.position:
            (lon, lat) = event.position
            posn_str = ('%.*f / %.*f'
                        % (LonLatPrecision, lon, LonLatPrecision, lat))
            fast_picture, slow_picture = \
              self.pyslip.tiles.lon_lat_to_picture_fast_slow(lon,lat)

            posn_str = ("Picture:  slow=%.*f / fast=%.*f pixels."
                        % (LonLatPrecision, slow_picture,
                           LonLatPrecision, fast_picture))
            coords = self.pyslip.tiles.get_flex_pixel_coordinates(lon,lat)
            if (len(coords) >= 2):
                if len(coords) == 3:
                    readout = int(round(coords[2]))
                else:
                    readout = -1

                coords_str = ("slow=%.*f / fast=%.*f pixels"
                              % (LonLatPrecision, coords[0],
                                 LonLatPrecision, coords[1]))
                if (len(coords) == 2):
                    posn_str += " Readout: " + coords_str +"."
                elif (readout >= 0):
                    posn_str += " Readout %d: %s." % (readout, coords_str)

                possible_intensity = None
                fi = self.pyslip.tiles.raw_image
                detector = fi.get_detector()
                ifs = (int(coords[1]), int(coords[0])) # int fast slow
                isf = (int(coords[0]), int(coords[1])) # int slow fast
                raw_data = fi.get_raw_data()
                if not isinstance(raw_data, tuple):
                  raw_data = (raw_data,)
                if len(detector) > 1:
                  if readout >= 0:
                    if detector[readout].is_coord_valid(ifs):
                        possible_intensity = raw_data[readout][isf]
                else:
                  if detector[0].is_coord_valid(ifs):
                      possible_intensity = raw_data[0][isf]

                if possible_intensity is not None:
                    if possible_intensity == 0:
                      format_str = " I=%6.4f"
                    else:
                      import math
                      yaya = int(math.ceil(math.log10(abs(possible_intensity))))
                      format_str = " I=%%6.%df"%(max(0,5-yaya))
                    posn_str += format_str%possible_intensity

                if (len(coords) > 2 and readout >= 0): # indicates it's a tiled image in a valid region
                    reso = self.pyslip.tiles.get_resolution(coords[1], coords[0], readout)
                else:
                    reso = self.pyslip.tiles.get_resolution(coords[1], coords[0])

                if reso is not None:
                    posn_str += " Resolution: %.3f"%(reso)

            self.statusbar.SetStatusText(posn_str)
        else:
            self.statusbar.SetStatusText("Click and drag to pan; "+
        "middle-click and drag to plot intensity profile, right-click to zoom")
            #print "event with no position",event
        return
        self.mouse_position.SetValue(posn_str)

    def handle_level_change(self, event):
        """Handle a pySlip LEVEL event."""
        pass
        return
        self.map_level.SetLabel('%d' % event.level)

    ######
    # Handle adding/removing select handler functions.
    ######

    def add_select_handler(self, id, handler):
        """Add handler for select in layer 'id'."""

        self.demo_select_dispatch[id] = handler

    def del_select_handler(self, id):
        """Remove handler for select in layer 'id'."""

        del self.demo_select_dispatch[id]

################################################################################

def run_demo():

    # start wxPython app
    app = wx.App()
    app_frame = AppFrame(tile_dir =
       "/Users/nksauter/xtalsoft/pyslip/pyslip-read-only/tiles")
    app_frame.Show()

##    import wx.lib.inspection
##    wx.lib.inspection.InspectionTool().Show()

    app.MainLoop()
