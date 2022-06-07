
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import plan_sad_experiment
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry

def exercise():
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
  #print result.representative_values[:-1]
  assert approx_equal(result.representative_values[:-1],
   [2.5, 12, 10880.374304954881, 3.8438000679016113, 97.77777777777779, 0.009, 0.9467311684652722, 0.9467311684652722, 0.7396676207890701, 0.7396676207890701, 20.9016001033, 98.3151904665], eps=0.01)
  assert (95 < result.representative_values[-2] < 100)
  # Insulin S-SAD
  with open("tst_plan_sad_experiment.fa", "w") as f:
    f.write("""
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
  #print result.representative_values[:-1]
  assert approx_equal(result.representative_values[:-1],
  [1.5, 12, 17126.515109651198, 0.5562999844551086, 97.77777777777779, 0.009, 0.5095890245676317, 0.5095890245676317, 0.616249995680551, 0.616249995680551, 17.36855053251538, 89.82473404335383], eps=0.01)
#  [2, 12, 7225.2485618841, 0.5562999844551086, 97.77777777777779, 0.009, 0.5095890245676317, 0.5095890245676317, 0.616249995680551, 0.616249995680551, 13.5861823956, 82.6264245311 ], eps=0.01)
  assert (86 < result.representative_values[-2] < 92)
  # now with worse resolution
  args = [
    "seq_file=tst_plan_sad_experiment.fa",
    "atom_type=S",
    "resolution=3.0",
    "wavelength=1.54"
  ]
  result = plan_sad_experiment.run(args=args, out=null_out())
  #print result.representative_values[:-1]
  assert approx_equal(result.representative_values[:-1],
  [5, 12, 462.41590796058244, 0.5562999844551086, 97.77777777777779, 0.009, 0.5095890245676317, 0.5095890245676317, 0.616249995680551, 0.616249995680551, 3.55312435872, 16.8013009611], eps=0.01)

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

if (__name__ == "__main__"):
  exercise()
  print("OK")
