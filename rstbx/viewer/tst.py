from __future__ import absolute_import, division, print_function

def exercise_1():
  from rstbx.viewer import screen_params
  p = screen_params()
  p.set_image_size(2048,2048)
  p.set_screen_size(1024,768)
  p.set_thumbnail_size(256,256,8)
  p.set_detector_resolution(0.1024)
  p.set_zoom(0)
  assert (p.get_scale() == 0.375)
  assert (p.get_bitmap_params() == (0, 0, 2048, 2048))
  p.set_zoom(0.5)
  assert (p.get_scale() == 0.5)
  assert (p.get_bitmap_params() == (0, 256, 2048, 1536))
  p.set_zoom(1.0)
  assert (p.get_bitmap_params() == (512, 640, 1024, 768))
  p.translate_image(-800, -800)
  assert (p.get_bitmap_params() == (1024, 1280, 1024, 768))
  p.translate_image(-400, -400)
  assert (p.get_bitmap_params() == (1024, 1280, 1024, 768))
  assert (p.get_thumbnail_box() == (128, 160, 128, 96))
  p.center_view_from_thumbnail(128,128)
  assert (p.get_bitmap_params() == (512, 640, 1024, 768))

# TODO: use a real image
def exercise_2():
  pass

if (__name__ == "__main__"):
  exercise_1()
  print("OK")
