
"""
Standalone tool for analyzing output of phenix.emringer and rendering plots,
as described in:

  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS. (2015)
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Nature Methods, in press.
"""

from __future__ import absolute_import, division, print_function
import mmtbx.ringer.em_scoring
import argparse
import os.path
import sys

def run(args, out=sys.stdout):
    parser = argparse.ArgumentParser()
    parser.add_argument("files",nargs="*")
    parser.add_argument("-s", "--Sampling_Angle", dest="sampling_angle", help="Don't mess with this unless you've also made the corresponding change in ringer. By default it is 5, which is identical to the default in ringer.", nargs='?', default=5)
    parser.add_argument("-r", "--Residues", dest="residues")
    parser.add_argument("--gui", dest="show_gui", action="store_true",
      default=False)
    args = parser.parse_args(args)
    #if (not args.show_gui):
    try :
      import matplotlib
    except ImportError as e :
      print("WARNING: matplotlib not present, plotting disabled", file=out)
      matplotlib = None
      args.show_gui = False
    else :
      matplotlib.use("Agg")
    app = None
    for file_name in args.files :
      result = mmtbx.ringer.em_scoring.main(
        file_name=file_name,
        sampling_angle=args.sampling_angle,
        out=out,
        quiet=(matplotlib is None)).show_summary(out=out)
      file_name = os.path.basename(file_name)
      if (args.show_gui):
        import wxtbx.plots.emringer
        import wxtbx.app
        if (app is None):
          app = wxtbx.app.CCTBXApp(0)
        f1 = wxtbx.plots.emringer.peaks_plot_frame(
          parent=None,
          title="Histogramss for %s" % file_name)
        f1.SetResult(result)
        f1.Show()
        f2 = wxtbx.plots.emringer.threshold_plot_frame(
          parent=None,
          title="Statistics across all thresholds for %s" % file_name)
        f2.SetResult(result)
        f2.Show()
    if (args.show_gui):
      app.MainLoop()

if __name__ == "__main__":
  run(sys.argv[1:])
