# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
from six.moves import range

import wx

from collections import OrderedDict
_scores = OrderedDict()


class ScoreSettingsFrame(wx.MiniFrame):
  # XXX Is this scoring or rating?

  def __init__(self, *args, **kwds):
    super(ScoreSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = ScoreSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class ScoreSettingsPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    from wxtbx import bitmaps

    super(ScoreSettingsPanel, self).__init__(*args, **kwds)

    # Needed for communication with the root frame.
    self._root_frame = self.GetParent().GetParent()

    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    self._id_text = wx.NewId()
    text = wx.StaticText(self, id= self._id_text, label="")
    sizer.Add(text, flag=wx.ALIGN_CENTER)

    # The score buttons.  XXX Would really like this to be stars or
    # some such.  Could then zap the static text above.
    box = wx.BoxSizer(wx.HORIZONTAL)
    self._id_buttons = []
    for i in range(6):
      btn = wx.Button(self, label="%d" % i)
      box.Add(btn, flag=wx.ALIGN_CENTER, border=5)
      self.Bind(wx.EVT_BUTTON, self.OnScore, btn)
      self._id_buttons.append(btn.GetId())
    sizer.Add(box)

    # XXX Would like to have label and short help for the buttons.
    # Buttons in wxWidgets version 2.9.1 and later support text and
    # bitmap, see SetBitmap, SetBitmapLabel, SetBitmapDisabled, etc.
    box = wx.BoxSizer(wx.HORIZONTAL)

    self._id_previous = wx.NewId()
    btn = wx.BitmapButton(
      self,
      id=self._id_previous,
      bitmap=bitmaps.fetch_icon_bitmap('actions', '1leftarrow'),
      style=wx.BORDER_NONE)
    box.Add(btn, flag=wx.ALIGN_CENTER, border=5)
    self.Bind(wx.EVT_BUTTON, self.OnPrevious, btn)

    btn = wx.Button(self, label="Save scores")
    box.Add(btn, flag=wx.ALIGN_CENTER, border=5)
    self.Bind(wx.EVT_BUTTON, self.OnSave, btn)

    self._id_next = wx.NewId()
    btn = wx.BitmapButton(
      self,
      id=self._id_next,
      bitmap=bitmaps.fetch_icon_bitmap('actions', '1rightarrow'),
      style=wx.BORDER_NONE)
    box.Add(btn, flag=wx.ALIGN_CENTER, border=5)
    self.Bind(wx.EVT_BUTTON, self.OnNext, btn)

    sizer.Add(box, flag=wx.ALIGN_CENTER)

    # Register update events for the dynamic widgets.
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateNext, id=self._id_next)
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdatePrevious, id=self._id_previous)
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateText, id=self._id_text)

    for i in range(self._root_frame.image_chooser.GetCount()):
      _scores[self._root_frame.get_key(self._root_frame.image_chooser.GetClientData(i))] = None

  def OnNext(self, event):
    self._root_frame.OnNext(event)


  def OnPrevious(self, event):
    self._root_frame.OnPrevious(event)


  def OnScore(self, event):
    score = self._id_buttons.index(event.EventObject.GetId())
    key = self._root_frame.get_key(self._root_frame.image_chooser.GetClientData(self._root_frame.image_chooser.GetSelection()))
    _scores[key] = score
    self.OnNext(event)


  def OnSave(self, event):
    dialog = wx.FileDialog(
      self,
      defaultDir='',
      message="Save scores",
      style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
      wildcard="Text files (*.txt)|*.txt")
    if dialog.ShowModal() == wx.ID_OK:
      path = dialog.GetPath()
      if (path != ''):
        stream = open(path, "w")
        for (key, score) in _scores.iteritems():
          if score is None:
            print("%s None" % (key), file=stream)
          else:
            print("%s %d" % (key, score), file=stream)
        stream.close()
        print("Dumped scores to", path)


  def OnUpdateNext(self, event):
    root_next = self._root_frame.toolbar.FindById(wx.ID_FORWARD)
    event.Enable(root_next.IsEnabled())


  def OnUpdatePrevious(self, event):
    root_previous = self._root_frame.toolbar.FindById(wx.ID_BACKWARD)
    event.Enable(root_previous.IsEnabled())


  def OnUpdateText(self, event):
    key = self._root_frame.GetTitle()
    if key in _scores:
      event.SetText("Previous score: %d" % _scores[key])
    else:
      event.SetText("Not previously scored")

    # Call self.Layout() to recenter the updated text.
    self.Layout()
