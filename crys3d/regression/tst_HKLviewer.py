
import libtbx.load_env, os.path, time
from crys3d.hklviewer import hklview_frame


def exercise1():
  file_name = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/reflection_files/2caz.mtz",
  test=os.path.isfile)
  myHKLview = hklview_frame.HKLViewFrame(verbose=0)
  myHKLview.LoadReflectionsFile(file_name)
  # slice expanded sphere of reflections with clip plane and translate to
  # the 46st plane of reflections along the K axis
  philstr = """
  clip_plane {
    normal_vector = "K-axis"
    is_assoc_real_space_vector = True
    clip_width = 1.184
    hkldist = 46
    normal_vector_length_scale = -1
  }
  viewer {
    data_array {
      label = 'F,SIGF'
      datatype = 'Amplitude'
    }
    show_vector = "['K-axis', True]"
    fixorientation = *vector None
  }
  hkls {
    expand_to_p1 = True
    expand_anomalous = True
  }
  """
  myHKLview.update_from_philstr(philstr)
  time.sleep(2)
  myHKLview.__exit__()
  refls = myHKLview.viewer.visible_hkls[:]
  # check that only the following 108 reflections were visible
  assert set(refls) == set( [(-24, 46, 1), (-17, 46, 0), (-28, 46, 1), (-24, 46, 0), (-17, 46, 1),
                    (-28, 46, 0), (-24, 46, -2), (-29, 46, 0), (-22, 46, 2), (-18, 46, -1),
                    (-26, 46, 2), (-24, 46, 2), (-29, 46, 1), (-22, 46, 1), (-26, 46, 1),
                    (-22, 46, 0), (-19, 46, -1), (-23, 46, 0), (-26, 46, 0), (-27, 46, 0),
                    (-23, 46, 1), (-21, 46, 2), (-27, 46, 1), (-25, 46, 2), (-20, 46, 1),
                    (-23, 46, 2), (-25, 46, -1), (-25, 46, -2), (-20, 46, 0), (-21, 46, -1),
                    (-21, 46, -2), (-21, 46, 0), (-25, 46, 0), (-21, 46, 1), (-25, 46, 1),
                    (-20, 46, 2), (-20, 46, -1), (-27, 46, -1), (-20, 46, -2), (-23, 46, -1),
                    (-23, 46, -2), (-19, 46, 0), (-19, 46, 1), (-26, 46, -1), (-26, 46, -2),
                    (-22, 46, -1), (-29, 46, -1), (-22, 46, -2), (-18, 46, 1), (-17, 46, -1),
                    (-18, 46, 0), (-28, 46, -1), (-24, 46, -1)]
              )

def run():
  exercise1()
  print("OK")

if __name__ == '__main__':
  run()
