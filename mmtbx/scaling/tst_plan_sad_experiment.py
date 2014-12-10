
from __future__ import division
from mmtbx.command_line import plan_sad_experiment
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry

def exercise () :
  # Generic SeMet protein (actually Rv0577)
  args = [
    "resolution=2.2",
    "atom_type=Se",
    "residues=300",
    "wavelength=0.9794",
    "include_weak_anomalous_scattering=False",
    "sites=12",
  ]
  result = plan_sad_experiment.run(args=args, out=null_out()).show(null_out())
  assert approx_equal(result.representative_values[:-1],
  [2.5, 12, 10880.374304954881, 3.8438000679016113, 97.77777777777779, 0.009, 0.5986974983459628, 0.8193812288687715, 0.6349909069671305, 0.7117990624917082, 25.471076479884275, 98.64197530864197], eps=0.01)
  assert (95 < result.representative_values[-2] < 100)
  # Insulin S-SAD
  open("tst_plan_sad_experiment.fa", "w").write("""
>1ZNI:A|PDBID|CHAIN|SEQUENCE
GIVEQCCTSICSLYQLENYCN
>1ZNI:B|PDBID|CHAIN|SEQUENCE
FVNQHLCGSHLVEALYLVCGERGFFYTPKA
>1ZNI:C|PDBID|CHAIN|SEQUENCE
GIVEQCCTSICSLYQLENYCN
>1ZNI:D|PDBID|CHAIN|SEQUENCE
FVNQHLCGSHLVEALYLVCGERGFFYTPKA
""")
  args = [
    "seq_file=tst_plan_sad_experiment.fa",
    "atom_type=S",
    "resolution=1.2",
    "wavelength=1.54"
  ]
  result = plan_sad_experiment.run(args=args, out=null_out())
  assert approx_equal(result.representative_values[:-1],
  [2, 12, 7225.2485618841, 0.5562999844551086, 97.77777777777779, 0.009, 0.3397568876945892, 0.6171343551165246, 0.5685560850183102, 0.6552289404785244, 11.647835631770215, 81.65725114921844], eps=0.01)
  assert (78 < result.representative_values[-2] < 84)
  # now with worse resolution
  args = [
    "seq_file=tst_plan_sad_experiment.fa",
    "atom_type=S",
    "resolution=3.0",
    "wavelength=1.54"
  ]
  result = plan_sad_experiment.run(args=args, out=null_out())
  assert approx_equal(result.representative_values[:-1],
  [5, 12, 462.41590796058244, 0.5562999844551086, 97.77777777777779, 0.009, 0.4438664886201641, 0.6584188448691838, 0.6044807178965403, 0.6683135637602184, 5.616135906033907, 34.428663031164064], eps=0.01)

  # Error handling
  args = [
    "resolution=2.2",
    "atom_type=Se",
    "wavelength=0.9794",
    "sites=12",
  ]
  try :
    result = plan_sad_experiment.run(args=args, out=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected

if (__name__ == "__main__") :
  exercise()
  print "OK"
