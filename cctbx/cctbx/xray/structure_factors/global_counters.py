import sys

calls_from_scatterers_direct = 0
time_from_scatterers_direct = 0
calls_from_scatterers_fft = 0
time_from_scatterers_fft = 0
calls_gradients_direct = 0
time_gradients_direct = 0
calls_gradients_fft = 0
time_gradients_fft = 0

def show(out=None, prefix=""):
  if (out is None): out = sys.stdout
  print >> out, prefix + "from_scatterers_direct: %d calls, %.2f s" % (
    calls_from_scatterers_direct, time_from_scatterers_direct)
  print >> out, prefix + "from_scatterers_fft:    %d calls, %.2f s" % (
    calls_from_scatterers_fft, time_from_scatterers_fft)
  print >> out, prefix + "gradients_direct:       %d calls, %.2f s" % (
    calls_gradients_direct, time_gradients_direct)
  print >> out, prefix + "gradients_fft:          %d calls, %.2f s" % (
    calls_gradients_fft, time_gradients_fft)
