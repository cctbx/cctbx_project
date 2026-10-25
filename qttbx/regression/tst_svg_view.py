"""SVG viewer tests: an SVG is shown whole and fitted, zooming leaves the
fitted state, and a file that is not SVG is refused."""

import os
import shutil
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.svg_view import SvgView, SvgViewer


# A tall picture, as a tree of many steps is: 200 wide, 2000 high.
_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="200" height="2000">
<line x1="100" y1="40" x2="100" y2="1900" stroke="black" stroke-width="2"/>
<polygon points="100,1900 93,1880 107,1880" fill="black"/>
<circle cx="100" cy="40" r="20" fill="orange" stroke="black"/>
<text x="100" y="74" text-anchor="middle" font-size="12">first</text>
</svg>
"""


def _qapp():
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  return app


def _shown(widget, width, height):
  # the viewport has its size only once the widget is shown
  widget.resize(width, height)
  widget.show()
  _qapp().processEvents()
  return widget


def _svg_file(tmp, text=_SVG, name="tree.svg"):
  path = os.path.join(tmp, name)
  with open(path, "w") as f:
    f.write(text)
  return path


def exercise_loads_whole_and_fitted():
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    view = _shown(SvgView(), 400, 400)
    assert view.load(_svg_file(tmp))
    rect = view.sceneRect()
    assert (rect.width(), rect.height()) == (200, 2000)
    # fitted: the whole height is inside the view, so the scale is below
    # a fifth and no scrolling is needed to reach the far end
    assert 0 < view.scale_factor() <= 400 / 2000.0
    corner = view.mapFromScene(rect.bottomRight())
    assert corner.y() <= view.viewport().height()
  finally:
    shutil.rmtree(tmp)


def exercise_zoom_leaves_fit_and_fit_returns():
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    view = _shown(SvgView(), 400, 400)
    view.load(_svg_file(tmp))
    fitted = view.scale_factor()
    view.zoom_in()
    assert abs(view.scale_factor() - fitted * SvgView.ZOOM_STEP) < 1e-6
    view.zoom_out()
    assert abs(view.scale_factor() - fitted) < 1e-6
    view.full_size()
    assert abs(view.scale_factor() - 1.0) < 1e-6
    # zoomed: a resize keeps the scale; fitted again: a resize refits
    _shown(view, 500, 500)
    assert abs(view.scale_factor() - 1.0) < 1e-6
    view.fit()
    _shown(view, 800, 800)
    assert view.scale_factor() > fitted
  finally:
    shutil.rmtree(tmp)


def exercise_scale_is_bounded():
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    view = _shown(SvgView(), 400, 400)
    view.load(_svg_file(tmp))
    for _ in range(100):
      view.zoom_in()
    assert abs(view.scale_factor() - SvgView.MAXIMUM_SCALE) < 1e-6
    for _ in range(200):
      view.zoom_out()
    assert abs(view.scale_factor() - SvgView.MINIMUM_SCALE) < 1e-6
  finally:
    shutil.rmtree(tmp)


def exercise_refuses_what_is_not_svg():
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    view = SvgView()
    assert not view.load(os.path.join(tmp, "missing.svg"))
    assert not view.load(_svg_file(tmp, text="not an svg", name="bad.svg"))
    assert len(view.scene().items()) == 0
  finally:
    shutil.rmtree(tmp)


def exercise_viewer_draws_the_picture():
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    viewer = _shown(SvgViewer(), 300, 600)
    assert viewer.load(_svg_file(tmp))
    assert viewer.windowTitle() == "tree.svg"
    # fitted, so the whole picture is on screen; every pixel is looked at,
    # as the line is thinner than any sampling step
    image = viewer.view.grab().toImage()
    white = orange = black = False
    for x in range(image.width()):
      for y in range(image.height()):
        c = image.pixelColor(x, y)
        r, g, b = c.red(), c.green(), c.blue()
        white = white or (r > 240 and g > 240 and b > 240)
        orange = orange or (r > 200 and 100 < g < 200 and b < 60)
        black = black or (r < 60 and g < 60 and b < 60)
    assert white and orange and black, (white, orange, black)
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_loads_whole_and_fitted()
  exercise_zoom_leaves_fit_and_fit_returns()
  exercise_scale_is_bounded()
  exercise_refuses_what_is_not_svg()
  exercise_viewer_draws_the_picture()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
