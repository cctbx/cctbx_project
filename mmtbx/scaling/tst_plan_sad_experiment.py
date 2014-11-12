
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
    [2.2, 12, 15965.989, 3.844, 97.778, 0.009, 0.509, 0.998, 0.617,
     0.774, 31.303,], eps=0.01)
  assert (95 < result.representative_values[-1] < 97)
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
  assert (not result.missed_target_resolutions)
  assert approx_equal(result.representative_values[:-1],
    [1.2, 12, 33450.225, 0.556, 97.778, 0.009, 0.393, 0.923, 0.591, 0.746,
     36.939,], eps=0.01)
  assert (96 < result.representative_values[-1] < 98)
  # now with worse resolution
  args = [
    "seq_file=tst_plan_sad_experiment.fa",
    "atom_type=S",
    "resolution=3.0",
    "wavelength=1.54"
  ]
  result = plan_sad_experiment.run(args=args, out=null_out())

  assert (result.missed_target_resolutions)
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
