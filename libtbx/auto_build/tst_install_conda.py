"""install_conda.py picks the Python interpreter of a conda environment.
Every platform must use the regular interpreter: bin/python on Linux and
macOS, python.exe on Windows. macOS used to run through the python.app
bundle, which GUI toolkits no longer need."""
import io
import os

from libtbx.auto_build.install_conda import conda_manager

def exercise_get_conda_python():
  prefix = os.path.join(os.sep, "root", "conda_base")
  m = conda_manager(root_dir=os.sep, conda_env=prefix, check_file=False,
                    log=io.StringIO())
  for system in ("Darwin", "Linux"):
    m.system = system
    python = m.get_conda_python(conda_env=prefix, check_file=False)
    assert python == os.path.join(prefix, "bin", "python"), (system, python)
  m.system = "Windows"
  python = m.get_conda_python(conda_env=prefix, check_file=False)
  assert python == os.path.join(prefix, "python.exe"), python

if __name__ == "__main__":
  exercise_get_conda_python()
  print("OK")
