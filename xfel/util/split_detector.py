from __future__ import division
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model.detector import Detector
from dials.array_family import flex
import sys

"""
Utility script for debugging. Example of invocation to split the panels
of a detector into 3x3 subpanels:
libtbx.python `libtbx.find_in_repositories xfel`/util/split_detector.py combined.expt combined.refl 3 3
"""

def split_detector(expts, refls, n_fast, n_slow):
  """
  Utitlity function to split a dxtbx detector model into subpanels.
  Each panel in the detector is split into n_fast x n_slow subpanels
  and reflections are reset to match their new panels.

  Columns affected: panel, xyzcal.mm, xyzcal.px, xyzobs.mm.value,
  xyzobs.mm.variance, xyzobs.px.value, and xyzobs.px.variance.
  """

  if "shoebox" in refls:
    print("Splitting shoeboxes not implemented, deleting them")
    del refls["shoebox"]

  assert len(expts.detectors()) == 1
  detector = expts.detectors()[0]
  new_detector = Detector()
  new_panel_ids = flex.size_t(len(refls))
  refl_f, refl_s, _ = refls["xyzobs.px.value"].parts()

  for panel_id, panel in enumerate(detector):
    image_fast, image_slow = panel.get_image_size()
    fast_axis = panel.get_fast_axis()
    slow_axis = panel.get_slow_axis()

    counter_fast = 0
    for fast in range(n_fast):
      counter_slow = 0
      for slow in range(n_slow):
        fast_pix = image_fast // n_fast + (1 if fast < image_fast % n_fast else 0)
        slow_pix = image_slow // n_slow + (1 if slow < image_slow % n_slow else 0)

        new_panel = new_detector.add_panel()
        new_panel.set_name(panel.get_name() + "_%d"%len(new_detector))
        new_panel.set_image_size((fast_pix, slow_pix))
        new_panel.set_frame(fast_axis, slow_axis, panel.get_pixel_lab_coord((counter_fast, counter_slow)))
        new_panel.set_pixel_size(panel.get_pixel_size())
        new_panel.set_trusted_range(panel.get_trusted_range())
        new_panel.set_thickness(panel.get_thickness())
        new_panel.set_material(panel.get_material())
        new_panel.set_mu(panel.get_mu())
        new_panel.set_gain(panel.get_gain())
        new_panel.set_px_mm_strategy(panel.get_px_mm_strategy())

        sel = (refl_f >= counter_fast) & (refl_f < counter_fast + fast_pix) & \
              (refl_s >= counter_slow) & (refl_s < counter_slow + slow_pix) & \
              (refls["panel"] == panel_id)
        subset = refls.select(sel)

        new_panel_ids.set_selected(sel, flex.size_t(len(subset), len(new_detector)-1))

        def reset_column(column):
          f, s, z = column.parts()
          f -= counter_fast; s -= counter_slow
          f_mm, s_mm = new_panel.pixel_to_millimeter(flex.vec2_double(f,s)).parts()
          return flex.vec3_double(f,s,z), flex.vec3_double(f_mm,s_mm,z)

        reset_px, reset_mm = reset_column(subset["xyzobs.px.value"])
        refls["xyzobs.px.value"].set_selected(sel, reset_px)
        refls["xyzobs.mm.value"].set_selected(sel, reset_mm)

        reset_px, reset_mm = reset_column(subset["xyzcal.px"])
        refls["xyzcal.px"].set_selected(sel, reset_px)
        refls["xyzcal.mm"].set_selected(sel, reset_mm)

        counter_slow += slow_pix
      counter_fast += fast_pix

  refls["panel"] = new_panel_ids
  for expt in expts:
    expt.detector = new_detector

  return expts, refls

if __name__ == "__main__":
  expts = ExperimentList.from_file(sys.argv[1], check_format=False)
  refls = flex.reflection_table.from_file(sys.argv[2])
  if len(sys.argv) > 3:
    n_fast, n_slow = map(int, sys.argv[3:5])
  else:
    n_fast = n_slow = 8
  expts, refls = split_detector(expts, refls, n_fast, n_slow)
  expts.as_file("split_%dx%d.expt"%(n_fast, n_slow))
  refls.as_file("split_%dx%d.refl"%(n_fast, n_slow))
