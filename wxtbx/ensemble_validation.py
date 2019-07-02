
"""
Classes for display of MolProbity validation statistics for multi-model PDB
files, used in GUI for phenix.ensemble_refinement.
"""

from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.utils import Sorry
from wxtbx import app, path_dialogs, plots
from mmtbx.command_line import validation_summary
from wxGUI2 import AdvancedWidgets, Base
import wx
import sys
from six.moves import range
from six.moves import zip

class ensemble_validation_plot(plots.histogram):
  def show_plot(self,
      values,
      as_histogram=False,
      n_bins=10,
      reference_value=None,
      title=None):
    if (as_histogram):
      self.show_histogram(
        data=values,
        n_bins=n_bins,
        reference_value=reference_value,
        x_label="Model Score",
        y_label="Number of Models",
        title=title)
    else :
      x = list(range(1, len(values) + 1))
      self.figure.clear()
      p = self.figure.add_subplot(111)
      ax = p.plot(x, values, '-^', color=(0.0,0.5,1.0))
      if (reference_value is not None):
        p.axhline(reference_value, color='red')
      p.set_xlabel("Model Number")
      p.set_ylabel("Score")
      if (title is not None):
        p.set_title(title)
      self.canvas.draw()

# =============================================================================
class ensemble_validation_panel(wx.Panel):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Plot controls
    box1 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box1)
    txt1 = wx.StaticText(self, label="Show statistic:")
    style = wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL
    box1.Add(txt1, 0, style , 5)
    self.stats_menu = wx.Choice(self,
      choices=validation_summary.molprobity_stat_labels)
    self.Bind(wx.EVT_CHOICE, self.OnSelectPlot, self.stats_menu)
    box1.Add(self.stats_menu, 0, style, 5)
    self.stats_menu.SetSelection(0)
    self.hist_box = wx.CheckBox(self, label="Display as histogram")
    box1.Add(self.hist_box, 0, style, 5)
    self.Bind(wx.EVT_CHECKBOX, self.OnSelectPlot, self.hist_box)
    n_bins_txt = wx.StaticText(self, label="Number of bins")
    box1.AddSpacer(40)
    box1.Add(n_bins_txt, 0, style, 5)
    self.n_bins_ctrl = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER,
                                   name="n_bins", size=(60, -1))
    self.n_bins_ctrl.SetValue("10")
    box1.Add(self.n_bins_ctrl, 0, style, 5)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnSelectPlot, self.n_bins_ctrl)
    self.n_bins_ctrl.Disable()
    self.save_ctrl = wx.Button(parent=self, label='Export raw data')
    box1.AddSpacer(40)
    box1.Add(self.save_ctrl, 0, style, 5)
    self.Bind(wx.EVT_BUTTON, self.OnSave, self.save_ctrl)

    # Data display
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box2)
    min_label = Base.BoldText(self, 'Minimum:')
    self.min_value = Base.PlainText(self, '')
    max_label = Base.BoldText(self, 'Maximum:')
    self.max_value = Base.PlainText(self, '')
    mean_label = Base.BoldText(self, 'Mean:')
    self.mean_value = Base.PlainText(self, '')
    self.value_ctrls = (self.min_value, self.max_value, self.mean_value)
    style = wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL
    for label, ctrl in zip([min_label, max_label, mean_label],self.value_ctrls):
      box2.Add(label, 0, style, 5)
      box2.Add(ctrl, 0, style, 5)
      box2.AddSpacer(40)
    self.plot = ensemble_validation_plot(
      parent=self,
      transparent=False)
    sizer.Add(self.plot, 1, wx.EXPAND|wx.ALL, 0)
    self.ensemble = None

  def set_ensemble(self, ensemble):
    assert (type(ensemble).__name__ == 'ensemble')
    self.ensemble = ensemble

  def OnSave(self, event):
    filename = path_dialogs.manager().select_file(
      parent=self,
      message='Save statistics as a CSV file',
      wildcard='CSV files (*.csv)|*.csv',
      current_file='statistics.csv',
      save=True)
    if (filename is not None):
      f = open(filename, 'w')
      labels = validation_summary.molprobity_stat_labels
      f.write('Statistic, Model 1, Model 2, ...\n')
      for label in labels:
        i_label = validation_summary.molprobity_stat_labels.index(label)
        stat = self.ensemble.__slots__[i_label]
        values = getattr(self.ensemble, stat)
        if ( (len(values) != 0) and (values.count(None) != len(values)) ):
          f.write('%s, ' % label +
                  ', '.join([str(v) for v in values]) + '\n')
      f.close()

  def OnSelectPlot(self, event):

    # plot line graph or histogram
    n_bins = 10
    as_histogram = self.hist_box.GetValue()
    if (as_histogram):
      self.n_bins_ctrl.Enable()
      n_bins = self.n_bins_ctrl.GetValue()
      try:
        n_bins = int(n_bins)
      except ValueError:
        raise Sorry('Please enter an integer for the number of histogram bins.')
    else:
      self.n_bins_ctrl.Disable()

    # check that n_bins is reasonable
    min_n_bins = 1
    max_n_bins = 1000
    if (n_bins < min_n_bins):
      n_bins = min_n_bins
    if (n_bins > max_n_bins):
      n_bins = max_n_bins
    self.n_bins_ctrl.SetValue(str(n_bins))

    # get attribute name and values for data
    label = self.stats_menu.GetStringSelection()
    i_label = validation_summary.molprobity_stat_labels.index(label)
    stat = self.ensemble.__slots__[i_label]
    values = getattr(self.ensemble, stat)
    if (len(values) == 0) or (values.count(None) == len(values)):
      self.plot.figure.clear()
      [ ctrl.SetLabel("") for ctrl in self.value_ctrls ]
      return
    mean = sum(values) / len(values)
    dist = [ min(values), max(values), mean ]
    [ ctrl.SetLabel("%.3f" % x) for x, ctrl in zip(dist, self.value_ctrls) ]
    self.plot.show_plot(
      values=values,
      as_histogram=as_histogram,
      n_bins=n_bins,
      reference_value=mean,
      title=label)
    self.Refresh()

# =============================================================================
class ensemble_chi_panel(wx.Panel):
  '''
  Automatic sequence of function calls
  load_data -> UpdateTable -> UpdatePlot
  OnSelect -> UpdatePlot

  load_data - initially loads data and sets up data
  UpdateTable - call this after changing threshold
  UpdatePlot - call this after changing number of histogram bins or autoscale
  '''
  def __init__(self, *args, **kwds):

    wx.Panel.__init__(self, *args, **kwds)
    self.main_window = self.GetParent().main_window
    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # controls
    control_sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(control_sizer)

    self.add_text(control_sizer, 'Threshold:')
    self.threshold_control = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER,
                                         name='threshold', size=(60,-1))
    self.add_control(control_sizer, self.threshold_control)
    self.Bind(wx.EVT_TEXT_ENTER, self.UpdateTable, self.threshold_control)

    self.add_text(control_sizer, 'Number of bins:')
    self.histogram_control = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER,
                                         name='n_bins', size=(60,-1))
    self.add_control(control_sizer, self.histogram_control)
    self.Bind(wx.EVT_TEXT_ENTER, self.UpdatePlot, self.histogram_control)

    self.autoscale_control = wx.CheckBox(self, label='Autoscale x-axis')
    self.add_control(control_sizer, self.autoscale_control)
    self.Bind(wx.EVT_CHECKBOX, self.UpdatePlot, self.autoscale_control)

    self.save_control = wx.Button(parent=self, label='Export raw data')
    self.add_control(control_sizer, self.save_control)
    self.Bind(wx.EVT_BUTTON, self.OnSave, self.save_control)

    # table
    self.table = AdvancedWidgets.OutlierList(
      parent=self,
      column_labels=[ 'Residue', 'Chi1', 'Chi2', 'Chi3', 'Chi4', 'Chi5' ],
      size=(600,100),
      column_formats=[ '%s', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f' ],
      style=wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.LC_VIRTUAL|wx.SUNKEN_BORDER)
    self.table.SetColumnWidths([200] + [100] * 5)
    sizer.Add(self.table, 0, wx.EXPAND)
    self.table.Bind(wx.EVT_LEFT_DOWN, self.OnSelect)

    # plot
    self.plot = plots.histogram(parent=self, transparent=False)
    sizer.Add(self.plot, 1, wx.EXPAND|wx.ALL, 0)

    # set defaults
    self.threshold_control.SetValue('0.95')
    self.histogram_control.SetValue('10')
    self.autoscale_control.SetValue(False)
    self.chi_angles = None
    self.meets_threshold = None
    self.row = 0
    self.col = 1

  def load_data(self, chi_angles=None):
    # reorganize data from
    # a list of all dihedrals for each atom_group per model
    # to
    # atom_group -> list of grouped dihedral angles (e.g. all chi1 in one list)
    if (chi_angles is not None):
      n_models = len(chi_angles)
      n_groups = len(chi_angles[0]['id_str'])

      id_str = chi_angles[0]['id_str']
      xyz = chi_angles[0]['xyz']
      id_str_map = dict()
      all_angles = [ list() for i in range(n_groups) ]
      for i in range(n_groups):           # loop over atom_groups
        n_dihedrals = len(chi_angles[0]['chi_angles'][i])
        angles = [ list() for j in range(n_dihedrals) ]
        for j in range(n_models):
          dihedrals = chi_angles[j]['chi_angles'][i]
          for k in range(n_dihedrals):
            dihedral = dihedrals[k]
            if (dihedral is not None):
              angles[k].append(dihedral)

        # check if adding 360 to negative angles is helpful
        # avoids issues where angles are clustered only near 0 and 360
        for j in range(len(angles)):
          if (len(angles[j]) > 0):
            stddev_a = flex.mean_and_variance(
              flex.double(angles[j])).unweighted_sample_variance()
            dihedrals_b = flex.double(angles[j])
            for k in range(len(dihedrals_b)):
              if (dihedrals_b[k] < 0.0):
                dihedrals_b[k] += 360
            stddev_b = flex.mean_and_variance(
              dihedrals_b).unweighted_sample_variance()
            if ( (stddev_b < stddev_a) or (max(angles[j]) < 0.0) ):
              angles[j] = list(dihedrals_b)

        # store reorganized data
        all_angles[i] = angles

        # build map matching id_str with row index for fast random access
        id_str_map[id_str[i]] = i

      self.chi_angles = { 'id_str': id_str,
                          'id_str_map': id_str_map,
                          'xyz': xyz,
                          'values': all_angles }
      self.meets_threshold = [ False for i in
                               range(len(self.chi_angles['id_str'])) ]
      self.UpdateTable()
      self.UpdatePlot()

  def OnSelect(self, event=None):
    self.row, self.col = self.table.GetRowCol(event)
    if (self.col == 0):
      self.col = 1
    self.UpdatePlot()

  def OnSave(self, event=None):
    filename = path_dialogs.manager().select_file(
      parent=self,
      message='Save dihedral angles as a CSV file',
      wildcard='CSV files (*.csv)|*.csv',
      current_file='dihedrals.csv',
      save=True)
    if (filename is not None):
      f = open(filename, 'w')
      id_str = self.chi_angles['id_str']
      dihedrals = self.chi_angles['values']
      for i in range(5):
        f.write('Chi %i\n' % (i+1))
        f.write('Residue, Model 1, Model 2, ...\n')
        for j in range(len(id_str)):
          if (i < len(dihedrals[j])):
            f.write('%s, ' % id_str[j] +
                    ', '.join([str(d) for d in dihedrals[j][i]]) + '\n')
        f.write('\n\n')
      f.close()

  def UpdateTable(self, event=None):
    '''
    Construct table of residues that satisfy threshold
    '''

    # check threshold
    threshold = self.threshold_control.GetValue()
    threshold_error = 'Please enter a fraction, between 0.0 and 1.0, for the threshold.'
    try:
      threshold = float(threshold)
    except ValueError:
      raise Sorry(threshold_error)
    if ( (threshold < 0.0) or (threshold > 1.0) ):
      raise Sorry(threshold_error)
    threshold = 1.0 - threshold

    # check if dihedrals lie within 2 standard deviations of the mean.
    # for Gaussian distributions, 95% of the sample is within 2 standard
    # deviations of the mean (default threshold value)
    # set atom_group to be displayed if fraction is below threshold.
    # actual check is n_outliers/n_total > (1 - input_threshold)
    for i in range(len(self.chi_angles['id_str'])): # loop over model
      self.meets_threshold[i] = False
      for j in range(len(self.chi_angles['values'][i])): # loop over group
        dihedrals = flex.double(self.chi_angles['values'][i][j])
        if (len(dihedrals) > 0):
          mean_stddev = flex.mean_and_variance(dihedrals)
          mean = mean_stddev.mean()
          stddev = mean_stddev.unweighted_sample_standard_deviation()
          n_outliers = 0
          chi_min = mean - 2.0*stddev
          chi_max = mean + 2.0*stddev
          for k in range(len(dihedrals)):
            if ( (dihedrals[k] < chi_min) or (dihedrals[k] > chi_max) ):
              n_outliers += 1
          is_outlier = (float(n_outliers)/len(dihedrals) > threshold)
          if (is_outlier):
            self.meets_threshold[i] = True
            break

    # refresh table
    table_data = list()
    for i in range(len(self.meets_threshold)):
      if (self.meets_threshold[i]):
        row = [ '---' for j in range(8) ]
        row[0] = self.chi_angles['id_str'][i]
        row[6] = None                          # sel_str for model viewer
        row[7] = self.chi_angles['xyz'][i]     # xyz for model viewer
        for j in range(len(self.chi_angles['values'][i])):
          chi = self.chi_angles['values'][i][j]
          if ( (None not in chi) and (len(chi) > 0) ):
            row[j+1] = flex.mean(flex.double(chi))
        table_data.append(row)
    self.table.ReloadData(table_data)
    self.table.Refresh()

  def UpdatePlot(self, event=None):
    '''
    Plot histogram of chi angle for residue and angle selected in table
    '''

    # check n_bins
    n_bins = self.histogram_control.GetValue()
    try:
      n_bins = int(n_bins)
    except ValueError:
      raise Sorry('Please enter an integer for the number of histogram bins.')

    # check that n_bins is reasonable
    min_n_bins = 1
    max_n_bins = 1000
    if (n_bins < min_n_bins):
      n_bins = min_n_bins
    if (n_bins > max_n_bins):
      n_bins = max_n_bins
    self.histogram_control.SetValue(str(n_bins))

    # check autoscale option
    autoscale = self.autoscale_control.GetValue()

    # update plot with data from currently selected cell
    if (True in self.meets_threshold):
      reference_value = self.table.get_item_text(self.row, self.col)
      if (reference_value != '---'):
        id_str = self.table.get_item_text(self.row, 0)
        actual_row = self.chi_angles['id_str_map'].get(id_str, 0)
        x_label = 'Chi %i (degrees)' % self.col
        y_label = 'Numer of Models'
        title = '%s Chi %i histogram' % (id_str, self.col)
        data=self.chi_angles['values'][actual_row][self.col-1]
        x_lim = (0.0, 360.0)
        if (autoscale):
          x_lim = None
        elif(min(data) < 0.0):
          x_lim = (-180.0, 180.0)
        self.plot.show_histogram(
          data=data, n_bins=n_bins, reference_value=float(reference_value),
          x_label=x_label, y_label=y_label, title=title, x_lim=x_lim)

  # convenience functions for placing text and widgets
  def add_text(self, sizer, label,
               style=wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL):
    text = wx.StaticText(self, label=label)
    sizer.Add(text, 0, style, 5)

  def add_control(self, sizer, control,
                  style=wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL,
                  spacer=40):
    sizer.Add(control, 0, style, 5)
    sizer.AddSpacer(spacer)

# =============================================================================
if (__name__ == "__main__"):
  result = validation_summary.run(sys.argv[1:])
  if (type(result).__name__ != 'ensemble'):
    raise Sorry("Not an ensemble, graphics not available.")
  app = app.CCTBXApp(0)
  frame = wx.Frame(None, -1, "Ensemble validation")
  szr = wx.BoxSizer(wx.VERTICAL)
  panel = ensemble_validation_panel(frame)
  panel.set_ensemble(result)
  szr.Add(panel, 1, wx.EXPAND|wx.ALL)
  frame.SetSizer(szr)
  szr.Fit(panel)
  frame.Fit()
  panel.OnSelectPlot(None)
  frame.Show()
  app.MainLoop()
