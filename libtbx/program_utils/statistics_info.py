from __future__ import absolute_import, division, print_function

# Rules of thumb for interpreting program statistics (used in Phenix GUI)

# format:
#   label min max quality strict_cutoff
threshold_values = [
  ("R-work", 0.4, 0.45, "fair", True),
  ("R-work", 0.55, 1, "bad", True),
  ("R-work", 0.0, 0.4, "good", True),
  ("R-free", 0.4, 0.45, "fair", True),
  ("R-free", 0.55, 1, "bad", True),
  ("R-free", 0.0, 0.4, "good", True),
  ("Skew", 0.15, 1, "good", True),
  ("RMS(bonds)", 0.0, 0.016, "good", True),
  ("RMS(bonds)", 0.016, 0.02, "fair", False),
  ("RMS(bonds)", 0.02, 999, "poor", False),
  ("RMS(angles)", 0.0, 1.6, "good", True),
  ("RMS(angles)", 1.6, 2.0, "fair", False),
  ("RMS(angles)", 2.0, 180, "poor", False),
  ("Ramachandran outliers", 0.0, 0.2, "good", False),
  ("Ramachandran outliers", 0.2, 0.5, "fair", False),
  ("Ramachandran outliers", 0.5, 1.0, "poor", False),
  ("Ramachandran outliers", 1.0, 100.0, "bad", False),
  ("Ramachandran favored", 98.0, 100.0, "good", False),
  ("Ramachandran favored", 95.0, 97.999, "fair", False),
  ("Ramachandran favored", 90.0, 94.999, "poor", False),
  ("Ramachandran favored", 0.0, 89.9999, "bad", False),
  ("Clashscore", 0.0, 19.999, "good", False),
  ("Clashscore", 20.0, 39.999, "fair", False),
  ("Clashscore", 40.0, 59.999, "poor", False),
  ("Clashscore", 60.0, 999, "bad", False),
  ("Rotamer outliers", 0.0, 1.999, "good", False),
  ("Rotamer outliers", 2.0, 5.0, "fair", False),
  ("Rotamer outliers", 5.0, 10.0, "poor", False),
  ("Rotamer outliers", 10.0, 100.0, "bad", False),
  # XXX backwards compatibility
  ("RMSbonds", 0.0, 0.016, "good", True),
  ("RMSbonds", 0.016, 0.02, "fair", False),
  ("RMSbonds", 0.02, 999, "poor", False),
  ("RMSangles", 0.0, 1.6, "good", True),
  ("RMSangles", 1.6, 2.0, "fair", False),
  ("RMSangles", 2.0, 180, "poor", False),
  ("CC", 0.0, 0.3, "bad", False),
  ("CC", 0.3, 0.5, "poor", False),
  ("CC", 0.5, 0.7, "fair", False),
  ("CC", 0.7, 1.0, "good", False),
  ("Model-map CC", 0.0, 0.3, "bad", False),
  ("Model-map CC", 0.3, 0.5, "poor", False),
  ("Model-map CC", 0.5, 0.7, "fair", False),
  ("Model-map CC", 0.7, 1.0, "good", False),
  ("Residues", 0, 0, "bad", False),
]

precisions = [
  ("R-work", 4),
  ("R-free", 4),
  ("RMS(bonds)", 4),
  ("RMS(angles)", 3),
  ("FOM", 3),
  ("Figure of Merit", 3),
  ("Skew", 2),
  ("Skewness", 2),
  ("Avg. B-factor", 2),
  ("Ramachandran favored", 2),
  ("Ramachandran outliers", 2),
  ("Rotamer outliers", 2),
  ("Clashscore", 1),
  # XXX backwards compatibility
  ("RMSbonds", 4),
  ("RMSangles", 3),
  ("CC", 2),
  ("Model-map CC", 2),
]

keys_and_labels = {
  "r_free" : "R-free",
  "r_work" : "R-work",
  "rms_bonds" : "RMS(bonds)",
  "rms_angles" : "RMS(angles)",
  "map_cc" : "Model-map CC",
  "cc" : "CC",
  "n_res" : "Residues",
  "rama_out" : "Ramachandran outliers",
  "rama_fav" : "Ramachandran favored",
  "rota_out" : "Rotamer outliers",
  "clashscore" : "Clashscore",
  "skew" : "Skew",
  "fom" : "Figure of merit",
  "rna_puckers" : "RNA pucker outliers",
  "rna_suites" : "RNA suite outliers",
  "rna_bonds" : "RNA bonds outliers",
  "rna_angles" : "RNA angles outliers",
}

precisions_dict_ = None
def get_format(stat_name, default="%g"):
  global precisions_dict_
  if (precisions_dict_ is None):
    precisions_dict_ = dict(precisions)
  if (stat_name in precisions_dict_):
    return "%%.%df" % precisions_dict_[stat_name]
  elif (stat_name in keys_and_labels):
    stat_label = keys_and_labels[stat_name]
    if (stat_label in precisions_dict_):
      return "%%.%df" % precisions_dict_[stat_label]
  return default
