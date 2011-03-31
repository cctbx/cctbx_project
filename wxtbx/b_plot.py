
from wxtbx import plots
import wx
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys

master_phil = libtbx.phil.parse("""
b_plot
  .short_caption = B-factor plot
  .caption = This utility will plot the B-factors for each protein or nucleic \
    acid chain in the input PDB file.  Waters and ligands will be ignored.  \
    A more comprehensive plot is available as part of the PHENIX validation \
    program, which incorporates real-space correlation to electron density, \
    adds markers for geometry outliers, and interfaces with Coot and PyMOL.
  .style = box auto_align caption_img:icons/crystal_project/32x32/mimetypes/spreadsheet.png
{
  pdb_file = None
    .type = path
    .short_caption = PDB file
    .style = bold file_type:pdb
  average_b_over = *residue mainchain sidechain
    .type = choice(multi=False)
    .short_caption = Average B-factors over
    .style = bold
  plot_range = *by_chain each_100_residues
    .type = choice(multi=False)
    .short_caption = Range of plot
    .caption = One_chain_at_a_time Every_100_residues
    .style = hidden
}
""")

class residue_info (object) :
  def __init__ (self,
                chain_id,
                resseq,
                icode,
                has_altconf,
                has_partocc,
                avg_b) :
    adopt_init_args(self, locals())

class analyze (object) :
  def __init__ (self, pdb_hierarchy, xray_structure, params, out=sys.stdout) :
    from cctbx import adptbx
    from scitbx.array_family import flex
    self.plot_range = params.plot_range
    self.chains = []
    self.residues = []
    b_isos = xray_structure.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.0)
    occ = pdb_hierarchy.atoms().extract_occ()
    model = pdb_hierarchy.models()[0]
    for chain in model.chains() :
      main_conf = chain.conformers()[0]
      is_na = main_conf.is_na()
      is_protein = main_conf.is_protein()
      if (not is_protein) and (not is_na) :
        print >> out, "Skipping chain '%s' - not protein or DNA/RNA." %chain.id
        continue
      self.chains.append(chain.id)
      self.residues.append([])
      for residue_group in chain.residue_groups() :
        n_conformers = len(residue_group.atom_groups())
        rg_i_seqs = residue_group.atoms().extract_i_seq()
        rg_occ = residue_group.atoms().extract_occ()
        if (params.average_b_over == "residue") :
          use_i_seqs = rg_i_seqs
        elif (params.average_b_over == "mainchain") :
          use_i_seqs = []
          if (is_protein) :
            for j_seq, atom in enumerate(residue_group.atoms()) :
              #alab = atom.fetch_labels()
              if (atom.name in [" N  ", " C  ", " CA ", " O  "]) :
                use_i_seqs.append(rg_i_seqs[j_seq])
          else :
            raise Sorry("Mainchain-only mode not supported for nucleic acids.")
        else :
          use_i_seqs = []
          if (is_protein) :
            for j_seq, atom in enumerate(residue_group.atoms()) :
              if (not atom.name in [" N  ", " C  ", " CA ", " O  "]) :
                use_i_seqs.append(rg_i_seqs[j_seq])
        if (len(use_i_seqs) > 0) :
          has_partocc = ((flex.min(occ.select(use_i_seqs)) < 1.0) and
                         (n_conformers == 1))
          res_info = residue_info(
            chain_id=chain.id,
            resseq=residue_group.resseq_as_int(),
            icode=residue_group.icode,
            has_altconf=(n_conformers > 1),
            has_partocc=has_partocc,
            avg_b=flex.mean(b_isos.select(use_i_seqs)))
          self.residues[-1].append(res_info)

  def make_plots (self, plot_range=None) :
    if (plot_range is None) :
      plot_range = self.plot_range
    import numpy
    plots = []
    altconf_val = 2 #max(min([ resi.avg_b for resi in self.residues ]) - 2, 0)
    if (plot_range == "by_chain") :
      for chain, residues in zip(self.chains, self.residues) :
        resid_start = ("%d%s" % (residues[0].resseq,residues[0].icode)).strip()
        resid_end = ("%d%s" % (residues[-1].resseq,residues[-1].icode)).strip()
        chain_vals = [] #numpy.array([])
        is_altconf = [] #numpy.array([])
        is_partocc = [] #numpy.array([])
        labels = []
        last_resseq = None
        for residue in residues :
          if (last_resseq is not None) :
            if (residue.resseq > (last_resseq + 1)) :
              gap_size = residue.resseq - last_resseq
              chain_vals.extend([numpy.NaN]* gap_size)
              is_altconf.extend([numpy.NaN] * gap_size)
              is_partocc.extend([numpy.NaN] * gap_size)
              labels.extend([None] * gap_size)
          if (residue.has_altconf) :
            is_altconf.append(altconf_val)
          else :
            is_altconf.append(numpy.NaN)
          if (residue.has_partocc) :
            is_partocc.append(altconf_val)
          else :
            is_partocc.append(numpy.NaN)
          chain_vals.append(residue.avg_b)
          labels.append(("%d%s" % (residue.resseq, residue.icode)).strip())
          last_resseq = residue.resseq
        chain_label = "Chain '%s' (%s - %s)" % (chain, resid_start, resid_end)
        plots.append((chain_label, chain_vals, is_altconf, is_partocc, labels))
    return plots

def run (args=(), params=None, out=sys.stdout) :
  pdb_file = params.b_plot.pdb_file
  from iotbx import file_reader
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs = pdb_in.file_object.xray_structure_simple()
  return analyze(pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    params=params.b_plot,
    out=out)

def show_plot_frame (result, parent=None) :
  frame = BPlotFrame(parent, -1, "B-factor plot")
  plots = result.make_plots()
  if (len(plots) == 0) :
    raise Sorry("No suitable chains found in PDB file.")
  frame.set_plot_data(plots)
  frame.Show()

class BPlotFrame (plots.plot_frame) :
  def draw_top_panel (self) :
    panel = wx.Panel(self)
    psizer = wx.BoxSizer(wx.VERTICAL)
    panel.SetSizer(psizer)
    txt1 = wx.StaticText(panel, -1, "Plot residues:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    self.chooser = wx.Choice(panel, -1, choices=[], size=(200,-1))
    self.Bind(wx.EVT_CHOICE, self.OnChoosePlot, self.chooser)
    bsizer = wx.BoxSizer(wx.HORIZONTAL)
    bsizer.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    bsizer.Add(self.chooser, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    psizer.Add(bsizer)
    self.top_panel = panel

  def create_plot_panel (self) :
    return b_plot_panel(
      parent=self,
      figure_size=(8,6))

  def set_plot_data (self, plots) :
    assert (len(plots) > 0)
    self._plot_data = plots
    all_b = []
    for avg_b in [ b for a,b,c,d,e in plots ] :
      all_b.extend(avg_b)
    self.plot_panel.set_limits(min(all_b) - 2, max(all_b) + 2)
    items = [ p[0] for p in plots ]
    self.chooser.SetItems(items)
    self.chooser.SetSelection(0)
    self.plot_panel.set_plot(self._plot_data[0])

  def OnChoosePlot (self, event) :
    assert hasattr(self, "_plot_data")
    choice = event.GetEventObject().GetSelection()
    self.plot_panel.set_plot(self._plot_data[choice])

def extract_labels (labels) :
  nx = len(labels)
  tickmarks = []
  show_labels = []
  if (nx <= 100) :
    for i, s in enumerate(labels) :
      if (s is not None) and (s.endswith("0")) :
        tickmarks.append(i+1)
        show_labels.append(s)
  elif (nx <= 250) :
    for i, s in enumerate(labels) :
      if (s is not None) and (s[-2:] in ["00", "25", "50", "75"]) :
        tickmarks.append(i+1)
        show_labels.append(s)
  elif (nx <= 500) :
    for i, s in enumerate(labels) :
      if (s is not None) and (s[-2:] in ["00", "50"]) :
        tickmarks.append(i+1)
        show_labels.append(s)
  else :
    for i, s in enumerate(labels) :
      if (s is not None) and (s.endswith("00")) :
        tickmarks.append(i+1)
        show_labels.append(s)
  return (tickmarks, show_labels)

class b_plot_panel (plots.plot_container) :
  def set_limits (self, ymin, ymax) :
    self._ymin = ymin
    self._ymax = ymax

  def set_plot (self, plot) :
    import numpy
    avg_b, is_alt, is_partial = plot[1], plot[2], plot[3]
    assert (len(avg_b) == len(is_alt) == len(is_partial))
    self.figure.clear()
    p = self.figure.add_subplot(1, 1, 1)
    a = p.get_axes()
    a.set_ylabel("B-factor")
    a.set_xlabel("Residue")
    p.set_title("Average B-factor: %s" % plot[0])
    x = range(1, len(plot[1]) + 1)
    p.plot(x, avg_b, "-", linewidth=1, color="r")
    xmarks, xlabels = extract_labels(plot[4])
    a.set_xticks(xmarks)
    a.set_xticklabels(xlabels)
    a.set_ylim(0, self._ymax)
    plot_labels = ["B(iso)"]
    if (set(is_alt) != set([numpy.NaN])) :
      p.plot(x, is_alt, "o")
      plot_labels.append("Alt. conf.")
    if (set(is_partial) != set([numpy.NaN])) :
      p.plot(x, is_partial, "^")
      plot_labels.append("Partial occupancy")
    self.figure.legend(p.lines, plot_labels, prop=self.get_font("legend"))
    self.canvas.draw()
    self.parent.Refresh()

def validate_params (params) :
  if (params.b_plot.pdb_file is None) :
    raise Sorry("No PDB file defined!")
