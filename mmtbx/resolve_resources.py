import libtbx.phil

phenix_masks_master_params = """\
solvent_content = 0.5
  .type = float
map_type = *2mFo-DFc
  .type = choice(multi=False)
resolution_factor = 0.25
  .type = float
probability_mask = True
  .type = bool
diff_map_cutoff = 1.5
  .type = float
output_all_masks = False
  .type = bool
"""

def get_phenix_masks_master_params():
  return libtbx.phil.parse(phenix_masks_master_params, process_includes=False)
