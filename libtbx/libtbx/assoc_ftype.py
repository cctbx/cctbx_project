import sys, os

def set():
  os.popen(r'assoc .px=PythonLibtbx', "r").read()
  py_exe = os.environ["LIBTBX_PYTHON_EXE"]
  os.popen(r'ftype PythonLibtbx="%s" "%%1" %%*' % py_exe, "r").read()

if (__name__ == "__main__"):
  set()
