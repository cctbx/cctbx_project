from __future__ import absolute_import, division, print_function
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex

# D. Rez, P. Rez & I. Grant
# Acta Cryst. (1994). A50, 481-497
table_2_stol = flex.double([
0.00,
0.05,
0.10,
0.15,
0.20,
0.25,
0.30,
0.35,
0.40,
0.45,
0.50,
0.60,
0.70,
0.80,
0.90,
1.00,
1.20,
1.40,
1.60,
1.80,
2.00,
2.50,
3.00,
3.50,
4.00,
5.00,
6.00])

table_2_o2minus = flex.double([
10.0000,
9.5884,
8.5345,
7.2188,
5.9558,
4.8934,
4.0574,
3.4190,
2.9364,
2.5720,
2.2963,
1.9275,
1.7085,
1.5678,
1.4644,
1.3774,
1.2191,
1.0668,
0.9215,
0.7881,
0.6696,
0.4404,
0.2907,
0.1951,
0.1332,
0.0686,
0.0351])

table_2_sigmas = flex.double(table_2_stol.size(), 0.00005)

def run():
  wk = xray_scattering.wk1995("O2-")
  gaussian_fit = scitbx.math.gaussian.fit(
    table_2_stol,
    table_2_o2minus,
    table_2_sigmas,
    wk.fetch())
  print("max error:", flex.max(gaussian_fit.significant_relative_errors()))
  cctbx.eltbx.gaussian_fit.write_plots(
    plots_dir="rez_plots",
    label=wk.label(),
    gaussian_fit=gaussian_fit)

if (__name__ == "__main__"):
  run()
