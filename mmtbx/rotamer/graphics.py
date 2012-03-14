
from libtbx.utils import Sorry
from libtbx import object_oriented_patterns as oop
from libtbx import adopt_init_args
import os

class rotarama_plot_mixin (object) :
  extent = [0, 360, 0, 360]
  def __init__ (self) :
    assert hasattr(self, "figure")
    self._points = []
    self._xyz = [] # only used by Phenix GUI (not offline plotting)
    self.plot = self.figure.add_subplot(111)
    self.plot.set_position([0.1, 0.1, 0.85, 0.85])

  def draw_plot (self,
                 stats,
                 title,
                 points=None,
                 show_labels=True,
                 colormap='jet',
                 contours=None,
                 xyz=None,
                 extent=None,
                 y_marks=None) :
    import matplotlib.cm
    self._points = []
    self._xyz = []
    cm = getattr(matplotlib.cm, colormap)
    self.plot.clear()
    if (extent is None) :
      extent = self.extent
    else :
      assert (len(extent) == 4)
    print extent
    self.plot.imshow(stats, origin="lower", cmap=cm, extent=extent)
    if (contours is not None) :
      self.plot.contour(stats, contours,
        origin="lower",
        colors='k',
        extent=extent)
    if (y_marks is None) :
      self.set_labels()
    else :
      self.set_labels(y_marks=y_marks)
    self.plot.set_title(title)
    if (points is not None) :
      if (xyz is not None) : assert (len(xyz) == len(points))
      for i, (x, y, label, is_outlier) in enumerate(points) :
        if is_outlier :
          self.plot.plot((x,),(y,), 'bo', markerfacecolor='red')
          if show_labels :
            self.plot.text(x, y, label, color='black')
          self._points.append((x,y))
          if (xyz is not None) :
            self._xyz.append(xyz[i])
        else :
          self.plot.plot((x,),(y,), 'bo', markerfacecolor='white')
    self.canvas.draw()

class ramachandran_plot_mixin (rotarama_plot_mixin) :
  extent = [-179,179,-179,179]
  def set_labels (self, y_marks=()) :
    axes = self.plot.get_axes()
    axes.set_xlabel("Phi")
    axes.set_xticks([-120,-60,0,60,120])
    axes.set_ylabel("Psi")
    axes.set_yticks([-120,-60,0,60,120])

class rotamer_plot_mixin (rotarama_plot_mixin) :
  def set_labels (self, y_marks=(60,180,300)) :
    axes = self.plot.get_axes()
    axes.set_xlabel("Chi1")
    axes.set_xticks([60,180,300])
    axes.set_ylabel("Chi2")
    axes.set_yticks(list(y_marks))
    axes.grid(True, color="0.75")

class simple_plot (object) :
  def __init__ (self,
                figure_size=(8,8),
                font_size=12,
                title_font_size=12,
                facecolor='white',
                transparent=False) :
    adopt_init_args(self, locals())
    try :
      import matplotlib
      import matplotlib.figure
      from matplotlib.backends.backend_agg import FigureCanvasAgg
    except ImportError, e :
      print e
      raise Sorry("Plotting requires that matplotlib be installed.")
    self.figure = matplotlib.figure.Figure(figure_size, 72, linewidth=0,
      facecolor=facecolor)
    if transparent :
      self.figure.figurePatch.set_alpha(0.0)
    self.canvas = FigureCanvasAgg(self.figure)
    #self.canvas.toolbar = oop.null()
    #self.figmgr = FigureManager(self.canvas, 1, self)

  def save_image (self, file_name) :
    assert (file_name is not None) and (file_name != "")
    base, ext = os.path.splitext(file_name)
    if (ext == ".pdf") :
      self.figure.savefig(file_name, orientation="landscape", format="pdf")
    elif (ext == ".ps") :
      self.figure.savefig(file_name, orientation="landscape", format="ps")
    elif (ext == ".png") :
      self.figure.savefig(file_name, format="png")
    else :
      raise RuntimeError("Extension %s not supported" % s)

class ramachandran_plot (simple_plot, ramachandran_plot_mixin) :
  def __init__ (self, *args, **kwds) :
    simple_plot.__init__(self, *args, **kwds)
    ramachandran_plot_mixin.__init__(self, *args, **kwds)

class rotamer_plot (simple_plot, rotamer_plot_mixin) :
  def __init__ (self, *args, **kwds) :
    simple_plot.__init__(self, *args, **kwds)
    rotamer_plot_mixin.__init__(self, *args, **kwds)

def get_residue_ramachandran_data (ramalyze_data,
                                   position_type,
                                   residue_name,
                                   point_type) :
  assert (position_type in ["general", "glycine", "cis-proline", "trans-proline",
    "pre-proline", "isoleucine or valine"])
  points, coords = [], []
  for i, residue in enumerate(ramalyze_data) :
    (chain_id,resseq,resname,quality,phi,psi,status,pos_name,xyz) = residue
    if (position_type == "general") :
      if (((residue_name == '*') or (resname.upper() == residue_name.upper()))
          and (not resname in ["PRO","GLY","ILE","VAL"])) :
        if ((point_type == "All") or
            (point_type=="Allowed/Outlier" and
             status in ["Allowed","OUTLIER"]) or
            (point_type == "Outlier" and status == "OUTLIER")) :
          points.append((phi, psi, "%s%s" % (chain_id, resseq),
            (status == "OUTLIER")))
          coords.append(xyz)
    elif (position_type.upper() == pos_name.upper()) :
      if ((point_type == "All") or
          (point_type=="Allowed/Outlier" and status in ["Allowed","OUTLIER"]) or
          (point_type == "Outlier" and status == "OUTLIER")) :
        points.append((phi, psi, "%s%s" % (chain_id, resseq),
          (status == "OUTLIER")))
        coords.append(xyz)
  return (points, coords)

def format_ramachandran_plot_title (position_type, residue_type) :
  title = "Ramachandran plot for "
  if (position_type == "cis-proline") :
    title += "cis-Proline"
  elif (position_type == "trans-proline") :
    title += "trans-Proline"
  elif (position_type == "glycine") :
    title += "Glycine"
  elif (position_type == "pre-proline") :
    title += "pre-Proline residues"
  elif (position_type == "isoleucine or valine") :
    title += "Ile or Val"
  else :
    if residue_type == '*' :
      title += "all non-Pro/Gly residues"
    else :
      title += residue_type
  return title

def draw_ramachandran_plot (ramalyze_data,
                            rotarama_data,
                            position_type,
                            file_name,
                            show_labels=True) :
  points, coords = get_residue_ramachandran_data(
    ramalyze_data=ramalyze_data,
    position_type=position_type,
    residue_name='*',
    point_type="All")
  p = ramachandran_plot()
  title = format_ramachandran_plot_title(position_type, '*')
  # XXX where do these numbers come from?
  if position_type == "general" :
    contours = [0.1495, 0.376]
  else :
    contours = [0.2115, 0.376]
  p.draw_plot(
    stats=rotarama_data,
    title=title,
    points=points,
    xyz=coords,
    colormap="Blues",
    contours=contours)
  p.save_image(file_name)

def get_residue_rotamer_data (rotalyze_data,
                              residue_name,
                              point_type) :
  points, coords = [], []
  for i, residue in enumerate(rotalyze_data) :
    (chain_id,resseq,resname,quality,chi1,chi2,chi3,chi4,status,xyz)=residue
    if resname.upper() == residue_name.upper() :
      if ((point_type == "All") or
          (point_type=="Outlier" and status=="OUTLIER")) :
        points.append((chi1, chi2, "%s%s" % (chain_id, resseq),
          (status == "OUTLIER")))
        coords.append(xyz)
  return (points, coords)

def draw_rotamer_plot (rotalyze_data,
                       rotarama_data,
                       residue_name,
                       file_name,
                       show_labels=True) :
  points, coords = get_residue_rotamer_data(
    rotalyze_data=rotalyze_data,
    residue_name=residue_name,
    point_type="All")
  p = rotamer_plot()
  title = "Chi1-Chi2 plot for %s" % residue_name
  p.draw_plot(
    stats=rotarama_data,
    title=title,
    points=points,
    xyz=coords,
    colormap="Blues",
    contours=None)
  p.save_image(file_name)
