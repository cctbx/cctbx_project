from __future__ import division
from mmtbx.refinement import tst_tardy_comprehensive
from mmtbx.refinement import tst_tardy_pdb
from scitbx.array_family import flex
from libtbx.math_utils import iceil
from libtbx import dict_with_default_0, group_args
from libtbx import easy_pickle
from libtbx import adopt_init_args
import sys, os
op = os.path

try: import reportlab
except ImportError: reportlab = None
if (reportlab is None):
  print "Skipping generation of plots: reportlab not available."

class plot_grid(object):

  def __init__(O,
        grid,
        top_labels,
        margin=20,
        top_label_space=80,
        more_narrow_shift=0):
    if (reportlab is None): return
    from reportlab.graphics.shapes import Group, String
    from reportlab.lib import pagesizes
    O.grid = grid
    O.margin = margin
    O.top_label_space = top_label_space
    O.more_narrow_shift = more_narrow_shift
    O.page_size = pagesizes.letter
    O.top_group = Group()
    for i,label in enumerate(top_labels):
      O.top_group.add(String(
        O.page_size[0] * 0.5,
        O.page_size[1] - O.margin - O.top_label_space * 0.3 * (i+1),
        label,
        fontSize=16,
        textAnchor="middle"))

  def process(O, grid_ij, xy_max, label, data, label_font_size=12):
    if (reportlab is None): return
    from reportlab.graphics.shapes import Group, String
    lp = O.line_plot(xy_max, data=data, label_font_size=label_font_size)
    gr = Group(lp)
    i,j = grid_ij
    assert 0 <= i < O.grid[0]
    assert 0 <= j < O.grid[1]
    i = O.grid[0] - 1 - i
    tx, ty = O.margin + j * (O.page_size[0] - 2 * O.margin) / O.grid[1] \
                      - j * O.more_narrow_shift, \
             O.margin + i * (O.page_size[1] - 2 * O.margin
                                            - O.top_label_space) / O.grid[0]
    gr.translate(tx, ty)
    O.top_group.add(gr)
    O.top_group.add(String(
      tx+lp.x+lp.width*0.5,
      ty+lp.y+lp.height*1.05,
      label,
      fontSize=label_font_size,
      textAnchor="middle"))
    if (i == 0 and j == 0):
      O.top_group.add(String(
        tx+lp.x+lp.width*0.5,
        ty+lp.y-lp.height*0.3,
        u"RMSD start (\u00C5)",
        fontSize=label_font_size,
        textAnchor="middle"))
      gr = Group(String(
        0,
        0,
        u"RMSD final (\u00C5)",
        fontSize=label_font_size,
        textAnchor="middle"))
      gr.rotate(90)
      gr.translate(
        ty+lp.y+lp.height*0.5,
        -(tx+lp.x-lp.width*0.15))
      O.top_group.add(gr)

  def line_plot(O, xy_max, data, label_font_size):
    if (reportlab is None): return
    from reportlab.graphics.charts.lineplots import LinePlot
    from reportlab.graphics.widgets.markers import makeMarker
    from reportlab.lib import colors
    lp = LinePlot()
    lp.x = 40
    lp.y = 40
    lp.height = 120
    lp.width = 120
    lp.data = [[(0,0),(xy_max,xy_max)], data]
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
    lp.xValueAxis.valueMax = xy_max
    lp.xValueAxis.valueSteps = range(xy_max+1)
    lp.xValueAxis.strokeWidth = 1
    lp.xValueAxis.tickDown = 3
    lp.xValueAxis.labels.fontSize = label_font_size
    lp.yValueAxis.valueMin = 0
    lp.yValueAxis.valueMax = xy_max
    lp.yValueAxis.valueSteps = range(xy_max+1)
    lp.yValueAxis.strokeWidth = 1
    lp.yValueAxis.tickLeft = 3
    lp.yValueAxis.labels.fontSize = label_font_size
    return lp

  def new_canvas(O, file_name):
    if (reportlab is None): return
    from reportlab.pdfgen.canvas import Canvas
    return Canvas(filename=file_name, pagesize=O.page_size)

  def draw_to_canvas(O, canvas):
    if (reportlab is None): return
    from reportlab.graphics.shapes import Drawing
    from reportlab.graphics import renderPDF
    drawing = Drawing(*O.page_size)
    drawing.add(O.top_group)
    renderPDF.draw(drawing=drawing, canvas=canvas, x=0, y=0)

  def write_to_file(O, file_name):
    if (reportlab is None): return
    canvas = O.new_canvas(file_name=file_name)
    O.draw_to_canvas(canvas=canvas)
    canvas.showPage()
    canvas.save()

class multi_page_plots(object):

  def __init__(O, file_name):
    if (reportlab is None): return
    O.file_name = file_name
    O.canvas = None

  def add_page(O, page):
    if (reportlab is None): return
    if (O.canvas is None):
      O.canvas = page.new_canvas(file_name=O.file_name)
    page.draw_to_canvas(canvas=O.canvas)
    O.canvas.showPage()

  def write_to_file(O):
    if (reportlab is None): return
    O.canvas.save()

def compose_top_label(pdb_file, random_displacements_parameterization, e):
  return ", ".join([
    pdb_file,
    "model: %s" % {False: "torsion", True: "cartesian"}[e],
    "random displacements: %s" % random_displacements_parameterization])

class min_mean_stats(object):

  def __init__(O, algorithm, random_displacements_parameterization, pdb_file):
    adopt_init_args(O, locals())
    O.data = {"tt": [], "tc": [], "ct": [], "cc": []}

  def collect(O, rmsd_t_c, param_values):
    mins = [flex.min(a) for a in rmsd_t_c]
    means = [flex.mean(a) for a in rmsd_t_c]
    if (mins[0] < mins[1]): a = "t"
    else:                   a = "c"
    if (means[0] < means[1]): b = "t"
    else:                     b = "c"
    O.data[a+b].append(param_values)

  def finalize(O):
    O.data = O.data.items()
    def cmp_data(a, b):
      result = -cmp(len(a[1]), len(b[1]))
      if (result == 0):
        result = cmp(a[0], b[0])
      return result
    O.data.sort(cmp_data)
    return O

  def show(O):
    for k,v in O.data:
      print "MIN_MEAN_STATS", O.algorithm, O.pdb_file, k, len(v), v

    print "MIN_MEAN_STATS", O.algorithm
    return O

  def pickle(O):
    easy_pickle.dump(file_name="min_mean_stats.pickle", obj=O)
    return O

def rmsd_start_final_plots_minimization(
      pdb_file,
      random_displacements_parameterization,
      tst_tardy_pdb_params,
      parameter_trial_table,
      cp_n_trials,
      rmsds,
      rmsd_n_n,
      write_separate_pages):
  ttd = dict(parameter_trial_table)
  plot_data = {}
  for h in ttd["structure_factors_high_resolution"]:
    plot_data[h] = {}
    for e in  ttd["emulate_cartesian"]:
      plot_data[h][e] = {}
      for d in ttd["real_space_gradients_delta_resolution_factor"]:
        plot_data[h][e][d] = {}
        for w in ttd["real_space_target_weight"]:
          plot_data[h][e][d][w] = group_args(
            rmsd_start=flex.double(),
            rmsd_final=flex.double(),
            rmsd_n=flex.double())
  #
  p = tst_tardy_pdb_params
  for cp_i_trial in xrange(cp_n_trials):
    tst_tardy_comprehensive.set_parameters(
      params=p,
      trial_table=parameter_trial_table,
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
  mpp = multi_page_plots(file_name="plots_h_e.pdf")
  for e in  ttd["emulate_cartesian"]:
    short_label = "e%d" % int(e)
    print short_label
    top_labels = [
      compose_top_label(
        pdb_file, random_displacements_parameterization, e),
      "algorithm: minimization"]
    page = plot_grid(grid=(4,3), top_labels=top_labels)
    w_d_ranks_rn = {}
    w_d_ranks_rf = {}
    for d in ttd["real_space_gradients_delta_resolution_factor"]:
      for w in ttd["real_space_target_weight"]:
        w_d_ranks_rn[(w,d)] = []
        w_d_ranks_rf[(w,d)] = []
    for i_h,h in enumerate(ttd["structure_factors_high_resolution"]):
      plot_xy_max = iceil(h * 1.3)
      page_rn = []
      page_rf = []
      for d in ttd["real_space_gradients_delta_resolution_factor"]:
        for i_w,w in enumerate(ttd["real_space_target_weight"]):
          pd = plot_data[h][e][d][w]
          rn = flex.mean(pd.rmsd_n)
          rf = flex.mean(pd.rmsd_final)
          page_rn.append((rn, (w,d)))
          page_rf.append((rf, (w,d)))
          label = "h%.2f_w%04.0f: %d=%.2f, %.2f" % (h, w, rmsd_n_n, rn, rf)
          print "  ", label
          page.process(
            grid_ij=(i_h,i_w),
            xy_max=plot_xy_max,
            label=label,
            data=zip(pd.rmsd_start, pd.rmsd_final))
      def cmp_rx(a, b):
        result = cmp(a[0], b[0])
        if (result == 0):
          result = cmp(a[1], b[1])
        return result
      page_rn.sort(cmp_rx)
      page_rf.sort(cmp_rx)
      for i,(r,w_d) in enumerate(page_rn):
        w_d_ranks_rn[w_d].append(i)
      for i,(r,w_d) in enumerate(page_rf):
        w_d_ranks_rf[w_d].append(i)
    if (write_separate_pages):
      page.write_to_file(file_name="plot_%s.pdf" % short_label)
    mpp.add_page(page=page)
    w_d_ranks_rn = w_d_ranks_rn.items()
    w_d_ranks_rf = w_d_ranks_rf.items()
    def cmp_w_d_ranks(a, b):
      result = cmp(sum(a[1]), sum(b[1]))
      if (result == 0):
        result = cmp(sorted(a[1]), sorted(b[1]))
        if (result == 0):
          result = cmp(a[1], b[1])
          if (result == 0):
            result = cmp(a[0], b[0])
      return result
    w_d_ranks_rn.sort(cmp_w_d_ranks)
    w_d_ranks_rf.sort(cmp_w_d_ranks)
    print "emulate_cartesian = %s" % str(e)
    for prefix,w_d_ranks in [("rn:", w_d_ranks_rn),
                             ("rf:", w_d_ranks_rf)]:
      for w_d,ranks in w_d_ranks:
        print prefix, "%4.0f %4.2f" % w_d, "%2d" % sum(ranks), \
          "[" + ", ".join(["%2d" % r for r in ranks]) + "]"
        prefix = "   "
      print
  mpp.write_to_file()
  #
  mms = min_mean_stats(
    algorithm="minimization",
    random_displacements_parameterization=random_displacements_parameterization,
    pdb_file=pdb_file)
  assert ttd["emulate_cartesian"] == (False, True)
  assert len(ttd["real_space_gradients_delta_resolution_factor"]) == 1
  for h in ttd["structure_factors_high_resolution"]:
    for d in ttd["real_space_gradients_delta_resolution_factor"]:
      for w in ttd["real_space_target_weight"]:
        mms.collect(
          rmsd_t_c=[plot_data[h][e][d][w].rmsd_final for e in (False, True)],
          param_values=(h,w))
  mms.finalize().show().pickle()

def rmsd_start_final_plots_annealing(
      pdb_file,
      random_displacements_parameterization,
      tst_tardy_pdb_params,
      parameter_trial_table,
      cp_n_trials,
      rmsds,
      write_separate_pages):
  ttd = dict(parameter_trial_table)
  plot_data = {}
  for h in ttd["structure_factors_high_resolution"]:
    plot_data[h] = {}
    for e in  ttd["emulate_cartesian"]:
      plot_data[h][e] = {}
      for d in ttd["real_space_gradients_delta_resolution_factor"]:
        plot_data[h][e][d] = {}
        for w in ttd["real_space_target_weight"]:
          plot_data[h][e][d][w] = {}
          for t in ttd["start_temperature_kelvin"]:
            plot_data[h][e][d][w][t] = {}
            for c in ttd["number_of_cooling_steps"]:
              plot_data[h][e][d][w][t][c] = group_args(
                rmsd_start=flex.double(),
                rmsd_final=flex.double())
  #
  p = tst_tardy_pdb_params
  for cp_i_trial in xrange(cp_n_trials):
    tst_tardy_comprehensive.set_parameters(
      params=p,
      trial_table=parameter_trial_table,
      cp_i_trial=cp_i_trial)
    plot = plot_data[
      p.structure_factors_high_resolution][
      p.emulate_cartesian][
      p.real_space_gradients_delta_resolution_factor][
      p.real_space_target_weight][
      p.start_temperature_kelvin][
      p.number_of_cooling_steps]
    for rmsd in rmsds[cp_i_trial]:
      plot.rmsd_start.append(rmsd[0])
      plot.rmsd_final.append(rmsd[-1])
  #
  extra_type = None
  if (pdb_file == "1yjp_box.pdb"):
    if (random_displacements_parameterization == "constrained"):
      extra_type = 0
  elif (pdb_file == "1yjp_no_water.pdb"):
    if (random_displacements_parameterization == "cartesian"):
      extra_type = 1
  if (extra_type is not None):
    extra_page = plot_grid(
      grid=(4,3), top_labels=[], more_narrow_shift=30)
  else:
    extra_page = None
  mpp = multi_page_plots(file_name="plots_h_e.pdf")
  for i_h,h in enumerate(ttd["structure_factors_high_resolution"]):
    for e in  ttd["emulate_cartesian"]:
      short_label = "h%.2f_e%d" % (h, int(e))
      print short_label
      top_labels = [compose_top_label(
        pdb_file, random_displacements_parameterization, e)]
      top_labels.append("high resol: %.2f, algorithm: annealing" % h)
      page = plot_grid(grid=(4,3), top_labels=top_labels)
      plot_xy_max = iceil(h * 1.3)
      for d in ttd["real_space_gradients_delta_resolution_factor"]:
        for i_w,w in enumerate(ttd["real_space_target_weight"]):
          for i_t,t in enumerate(ttd["start_temperature_kelvin"]):
            for i_c,c in enumerate(ttd["number_of_cooling_steps"]):
              pd = plot_data[h][e][d][w][t][c]
              rf = flex.mean(pd.rmsd_final)
              label = "w%04.0f_t%.0f_c%d: %.2f" % (w, t, c, rf)
              print "  ", label
              page.process(
                grid_ij=(i_t*2+i_c, i_w),
                xy_max=plot_xy_max,
                label=label,
                data=zip(pd.rmsd_start, pd.rmsd_final))
              if (extra_type == 0
                  and h == 3.75
                  and t == 5000
                  and c == 500):
                extra_label = "w_rs = %.0f" % w
                extra_page.process(
                  grid_ij=(i_w+1, int(e)),
                  xy_max=plot_xy_max,
                  label=extra_label,
                  data=zip(pd.rmsd_start, pd.rmsd_final),
                  label_font_size=14)
              elif (extra_type == 1
                    and w == 100
                    and (h == 3.75 or h == 5.00)
                    and c == 500):
                extra_label = u"resol. = %.2f \u00C5" % h
                extra_page.process(
                  grid_ij=(i_h, int(e)),
                  xy_max=plot_xy_max,
                  label=extra_label,
                  data=zip(pd.rmsd_start, pd.rmsd_final),
                  label_font_size=14)
      if (write_separate_pages):
        page.write_to_file(file_name="plot_%s.pdf" % short_label)
      mpp.add_page(page=page)
  mpp.write_to_file()
  if (extra_page is not None):
    from reportlab.graphics.shapes import String
    if (extra_type == 0): ty = 540
    else:                 ty = 372
    extra_page.top_group.add(String(
      120, ty,
      "Torsion-Angle SA",
      fontSize=16,
      textAnchor="middle"))
    extra_page.top_group.add(String(
      280+2/3, ty,
      "Cartesian SA",
      fontSize=16,
      textAnchor="middle"))
    extra_page.write_to_file(file_name="plot_extra.pdf")
  #
  mms = min_mean_stats(
    algorithm="annealing",
    random_displacements_parameterization=random_displacements_parameterization,
    pdb_file=pdb_file)
  assert ttd["emulate_cartesian"] == (False, True)
  assert len(ttd["real_space_gradients_delta_resolution_factor"]) == 1
  for h in ttd["structure_factors_high_resolution"]:
    for d in ttd["real_space_gradients_delta_resolution_factor"]:
      for w in ttd["real_space_target_weight"]:
        for t in ttd["start_temperature_kelvin"]:
          for c in ttd["number_of_cooling_steps"]:
            mms.collect(
              rmsd_t_c=[plot_data[h][e][d][w][t][c].rmsd_final
                for e in (False, True)],
              param_values=(h,w,t,c))
  mms.finalize().show().pickle()

def eval_1(args):
  first_file = open(args[0]).read().splitlines()
  #
  for line in first_file:
    if (line.startswith("pdb_file = ")):
      pdb_file = line.split(" ", 2)[2]
      assert pdb_file[0] == pdb_file[-1]
      assert pdb_file[0] in ['"', "'"]
      pdb_file = op.basename(pdb_file[1:-1])
      break
  else:
    raise RuntimeError('pdb_file = "..." not found.')
  #
  for line in first_file:
    if   (line ==
            "random_displacements_parameterization = *constrained cartesian"):
      random_displacements_parameterization = "constrained"
      break
    elif (line ==
            "random_displacements_parameterization = constrained *cartesian"):
      random_displacements_parameterization = "cartesian"
      break
  else:
    raise RuntimeError(
      "random_displacements_parameterization = constrained or cartesian"
      " not found.")
  #
  for line in first_file:
    if (line == "algorithm = *minimization annealing"):
      algorithm = "minimization"
      break
    elif (line == "algorithm = minimization *annealing"):
      algorithm = "annealing"
      break
  else:
    raise RuntimeError("algorithm = minimization or annealing not found.")
  #
  del first_file
  #
  tst_tardy_pdb_master_phil = tst_tardy_pdb.get_master_phil()
  tst_tardy_pdb_params = tst_tardy_pdb_master_phil.extract()
  if (algorithm == "minimization"):
    parameter_trial_table \
      = tst_tardy_comprehensive.common_parameter_trial_table
  elif (algorithm == "annealing"):
    parameter_trial_table \
      = tst_tardy_comprehensive.annealing_parameter_trial_table
  else:
    raise AssertionError
  cp_n_trials = tst_tardy_comprehensive.number_of_trials(
    table=parameter_trial_table)
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
  write_separate_pages = False
  if (algorithm == "minimization"):
    rmsd_start_final_plots_minimization(
      pdb_file=pdb_file,
      random_displacements_parameterization
        =random_displacements_parameterization,
      tst_tardy_pdb_params=tst_tardy_pdb_params,
      parameter_trial_table=parameter_trial_table,
      cp_n_trials=cp_n_trials,
      rmsds=rmsds,
      rmsd_n_n=50,
      write_separate_pages=write_separate_pages)
  elif (algorithm == "annealing"):
    rmsd_start_final_plots_annealing(
      pdb_file=pdb_file,
      random_displacements_parameterization
        =random_displacements_parameterization,
      tst_tardy_pdb_params=tst_tardy_pdb_params,
      parameter_trial_table=parameter_trial_table,
      cp_n_trials=cp_n_trials,
      rmsds=rmsds,
      write_separate_pages=write_separate_pages)

def eval_2(args):
  # HIGH REDUNDANCY eval_2_orca()
  i_row_by_tc = {
    "tt": 0,
    "tc": 1,
    "ct": 2,
    "cc": 3}
  i_col_major = {
    "gly_gly_box.pdb": 0,
    "lys_pro_trp_box.pdb": 1,
    "1yjp_box.pdb": 2,
    "1yjp_no_water.pdb": 3}
  i_col_minor = {
    "constrained": 0,
    "cartesian": 1}
  n_cols = len(i_col_major) * len(i_col_minor)
  i_table = {
    "minimization": 0,
    "annealing": 1}
  def make_table():
    return [[None] * n_cols for tc in i_row_by_tc]
  tabs = [make_table() for i_tab in xrange(4)]
  done = set()
  for file_name in args:
    mms = easy_pickle.load(file_name=file_name)
    i_tab = i_table[mms.algorithm]
    i_col = \
      i_col_major[mms.pdb_file] * len(i_col_minor) + \
      i_col_minor[mms.random_displacements_parameterization]
    for tc,list_of_param_values in mms.data:
      i_row = i_row_by_tc[tc]
      key = (i_tab, i_row, i_col)
      assert key not in done
      done.add(key)
      tabs[i_tab*2][i_row][i_col] = len(list_of_param_values)
      n = 0
      if (mms.algorithm == "minimization"):
        for h,w in list_of_param_values:
          if (w == 10): n += 1
      elif (mms.algorithm == "annealing"):
        for h,w,t,c in list_of_param_values:
          if (w == 10): n += 1
      else:
        raise AssertionError
      tabs[i_tab*2+1][i_row][i_col] = len(list_of_param_values) - n
  assert len(done) == len(tabs)//2 * len(i_row_by_tc) * n_cols
  table_i = dict([tuple(reversed(item)) for item in i_table.items()])
  for i_tab,tab in enumerate(tabs):
    print table_i[i_tab//2], "omit_weight_10=%s" % str(bool(i_tab%2))
    for row in tab:
      print " ".join([str(v) for v in row])
    print

def eval_2_orca(args):
  # HIGH REDUNDANCY eval_2()
  i_row_by_tc = {
    "tt": 0,
    "tc": 1,
    "ct": 2,
    "cc": 3}
  i_col_major = {
    "gly_gly_box.pdb": 0,
    "lys_pro_trp_box.pdb": 1,
    "1yjp_box.pdb": 2,
    "1yjp_no_water.pdb": 3}
  i_col_minor = {
    "constrained": 0}
  n_cols = len(i_col_major) * len(i_col_minor)
  i_table = {
    "minimization": 0}
  def make_table():
    return [[None] * n_cols for tc in i_row_by_tc]
  tabs = [make_table() for i_tab in xrange(2)]
  done = set()
  for file_name in args:
    mms = easy_pickle.load(file_name=file_name)
    i_tab = i_table[mms.algorithm]
    i_col = \
      i_col_major[mms.pdb_file] * len(i_col_minor) + \
      i_col_minor[mms.random_displacements_parameterization]
    for tc,list_of_param_values in mms.data:
      i_row = i_row_by_tc[tc]
      key = (i_tab, i_row, i_col)
      assert key not in done
      done.add(key)
      tabs[i_tab*2][i_row][i_col] = len(list_of_param_values)
      n = 0
      if (mms.algorithm == "minimization"):
        for h,w in list_of_param_values:
          if (w == 10): n += 1
      elif (mms.algorithm == "annealing"):
        for h,w,t,c in list_of_param_values:
          if (w == 10): n += 1
      else:
        raise AssertionError
      tabs[i_tab*1+1][i_row][i_col] = len(list_of_param_values) - n
  assert len(done) == len(tabs)//2 * len(i_row_by_tc) * n_cols
  table_i = dict([tuple(reversed(item)) for item in i_table.items()])
  for i_tab,tab in enumerate(tabs):
    print table_i[i_tab//2], "omit_weight_10=%s" % str(bool(i_tab%2))
    for row in tab:
      print " ".join([str(v) for v in row])
    print

def run(args):
  if (args[0] == "eval_pickles"):
    eval_2(args=args[1:])
  elif (args[0] == "eval_pickles_orca"):
    eval_2_orca(args=args[1:])
  else:
    eval_1(args=args)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
