from __future__ import division
from mmtbx.refinement import tst_tardy_comprehensive
from mmtbx.refinement import tst_tardy_pdb
from scitbx.array_family import flex
from libtbx import dict_with_default_0, group_args
import sys

try: import reportlab
except ImportError: reportlab = None

class plot_grid(object):

  def __init__(O, grid, top_label, margin=20, top_label_space=40):
    from reportlab.graphics.shapes import Group, String
    from reportlab.lib import pagesizes
    O.grid = grid
    O.margin = margin
    O.top_label_space = top_label_space
    O.page_size = pagesizes.letter
    O.top_group = Group()
    O.top_group.add(String(
      O.page_size[0] * 0.5,
      O.page_size[1] - O.margin - O.top_label_space * 0.6,
      top_label,
      fontSize=16,
      textAnchor="middle"))

  def process(O, grid_ij, label, data):
    from reportlab.graphics.shapes import Group, String
    lp = O.line_plot(data=data)
    gr = Group(lp)
    i,j = grid_ij
    tx, ty = O.margin + i * (O.page_size[0] - 2 * O.margin) / O.grid[0], \
             O.margin + j * (O.page_size[1] - 2 * O.margin
                                            - O.top_label_space) / O.grid[1]
    gr.translate(tx, ty)
    O.top_group.add(gr)
    f = flex.mean(flex.double(data.rmsd_final))
    n = flex.mean(flex.double(data.rmsd_n))
    info = "%s: %d=%.2f, %.2f" % (label, data.rmsd_n_n, n, f)
    O.top_group.add(String(
      tx+lp.x+lp.width*0.5,
      ty+lp.y+lp.height*1.05,
      info,
      fontSize=12,
      textAnchor="middle"))

  def line_plot(O, data):
    from reportlab.graphics.charts.lineplots import LinePlot
    from reportlab.graphics.widgets.markers import makeMarker
    from reportlab.lib import colors
    lp = LinePlot()
    lp.x = 40
    lp.y = 40
    lp.height = 120
    lp.width = 120
    lp.data = [[(0,0),(4,4)], zip(data.rmsd_start, data.rmsd_final)]
    lp.lines[0].strokeColor = colors.Color(*[0.8]*3)
    lp.lines[0].strokeWidth = 0.5
    lp.lines[1].strokeColor = colors.white
    lp.lines[1].strokeWidth = 0
    lp.lines[1].symbol = makeMarker("Circle")
    lp.lines[1].symbol.strokeWidth = 0.5
    lp.joinedLines = 1
    lp.strokeColor = colors.Color(*[0.8]*3)
    lp.strokeWidth = 1
    lp.xValueAxis.valueMin = 0
    lp.xValueAxis.valueMax = 4
    lp.xValueAxis.valueSteps = range(5)
    lp.xValueAxis.strokeWidth = 1
    lp.xValueAxis.tickDown = 3
    lp.xValueAxis.labels.fontSize = 10
    lp.yValueAxis.valueMin = 0
    lp.yValueAxis.valueMax = 4
    lp.yValueAxis.valueSteps = range(5)
    lp.yValueAxis.strokeWidth = 1
    lp.yValueAxis.tickLeft = 3
    lp.yValueAxis.labels.fontSize = 10
    return lp

  def new_canvas(O, file_name):
    from reportlab.pdfgen.canvas import Canvas
    return Canvas(filename=file_name, pagesize=O.page_size)

  def draw_to_canvas(O, canvas):
    from reportlab.graphics.shapes import Drawing
    from reportlab.graphics import renderPDF
    drawing = Drawing(*O.page_size)
    drawing.add(O.top_group)
    renderPDF.draw(drawing=drawing, canvas=canvas, x=0, y=0)

  def write_to_file(O, file_name):
    canvas = O.new_canvas(file_name=file_name)
    O.draw_to_canvas(canvas=canvas)
    canvas.showPage()
    canvas.save()

class multi_page_plots(object):

  def __init__(O, file_name):
    O.file_name = file_name
    O.canvas = None

  def add_page(O, page):
    if (O.canvas is None):
      O.canvas = page.new_canvas(file_name=O.file_name)
    page.draw_to_canvas(canvas=O.canvas)
    O.canvas.showPage()

  def write_to_file(O):
    O.canvas.save()

def rmsd_start_final_plots(
      tst_tardy_pdb_params,
      cp_n_trials,
      rmsds,
      rmsd_n_n,
      write_separate_pages):
  ttd = dict(tst_tardy_comprehensive.common_parameter_trial_table)
  plot_data = {}
  for h in ttd["structure_factors_high_resolution"]:
    plot_data[h] = {}
    for e in  ttd["emulate_cartesian"]:
      plot_data[h][e] = {}
      for d in ttd["real_space_gradients_delta_resolution_factor"]:
        plot_data[h][e][d] = {}
        for w in ttd["real_space_target_weight"]:
          plot_data[h][e][d][w] = group_args(
            rmsd_start=[],
            rmsd_final=[],
            rmsd_n=[],
            rmsd_n_n=rmsd_n_n)
  #
  p = tst_tardy_pdb_params
  for cp_i_trial in xrange(cp_n_trials):
    tst_tardy_comprehensive.set_parameters(
      params=p,
      trial_table=tst_tardy_comprehensive.common_parameter_trial_table,
      cp_i_trial=cp_i_trial)
    plot = plot_data[
      p.structure_factors_high_resolution][
      p.emulate_cartesian][
      p.real_space_gradients_delta_resolution_factor][
      p.real_space_target_weight]
    for rmsd in rmsds[cp_i_trial]:
      plot.rmsd_start.append(rmsd[0])
      plot.rmsd_final.append(rmsd[-1])
      plot.rmsd_n.append(rmsd[min(len(rmsd)-1, rmsd_n_n)])
  #
  if (reportlab is None):
    print "Skipping generation of plots: reportlab not available."
    return
  #
  mpp = multi_page_plots(file_name="plots_h_e.pdf")
  for h in ttd["structure_factors_high_resolution"]:
    for e in  ttd["emulate_cartesian"]:
      top_label = "h%.0f_e%d" % (h, int(e))
      page = plot_grid(grid=(3,4), top_label=top_label)
      for i,d in enumerate(ttd["real_space_gradients_delta_resolution_factor"]):
        for j,w in enumerate(ttd["real_space_target_weight"]):
          page.process(
            grid_ij=(i,j),
            label="w%04.0f_d%.0f" % (w, d*100),
            data=plot_data[h][e][d][w])
      if (write_separate_pages):
        page.write_to_file(file_name="plot_%s.pdf" % top_label)
      mpp.add_page(page=page)
  mpp.write_to_file()

def run(args):
  tst_tardy_pdb_master_phil = tst_tardy_pdb.get_master_phil()
  tst_tardy_pdb_params = tst_tardy_pdb_master_phil.extract()
  cp_n_trials = tst_tardy_comprehensive.number_of_trials(
    table=tst_tardy_comprehensive.common_parameter_trial_table)
  #
  random_seed_rmsd = []
  for cp_i_trial in xrange(cp_n_trials):
    random_seed_rmsd.append({})
  for file_name in args:
    for line in open(file_name).read().splitlines():
      if (not line.startswith("RESULT_cp_i_trial_random_seed_rmsd: ")):
        continue
      flds = line.split(None, 3)
      assert len(flds) == 4
      cp_i_trial = int(flds[1])
      random_seed = int(flds[2])
      rmsd = flex.double(eval(flds[3]))
      assert not random_seed_rmsd[cp_i_trial].has_key(random_seed)
      random_seed_rmsd[cp_i_trial][random_seed] = rmsd
  random_seeds_found = dict_with_default_0()
  for cp_i_trial in xrange(cp_n_trials):
    random_seeds_found[tuple(sorted(random_seed_rmsd[cp_i_trial].keys()))] += 1
  if (len(random_seeds_found) != 1):
    print random_seeds_found
    raise RuntimeError("Unexpected random_seeds_found (see output).")
  assert random_seeds_found.values()[0] == cp_n_trials
  random_seeds_found = random_seeds_found.keys()[0]
  assert random_seeds_found == tuple(range(len(random_seeds_found)))
  rmsds = []
  for cp_i_trial in xrange(cp_n_trials):
    rmsds.append([random_seed_rmsd[cp_i_trial][random_seed]
      for random_seed in xrange(len(random_seeds_found))])
  #
  if (1):
    rmsd_start_final_plots(
      tst_tardy_pdb_params=tst_tardy_pdb_params,
      cp_n_trials=cp_n_trials,
      rmsds=rmsds,
      rmsd_n_n=50,
      write_separate_pages=False)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
