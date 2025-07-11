'''
Data labels used in CIF and MTZ formats and a mapping
between the two formats

References:
http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/refln.html
http://www.ccp4.ac.uk/html/cif2mtz.html
'''
from __future__ import absolute_import, division, print_function

# =============================================================================
# intensities
phenix_to_cif_intensities = {
  # 'IOBS': '_refln.F_squared_meas',
  'IOBS': '_refln.intensity_meas',           # This is used in PDB
  # 'SIGIOBS': '_refln.F_squared_sigma',
  'SIGIOBS': '_refln.intensity_sigma',       # This is used in PDB
  'IOBS(+)': '_refln.pdbx_I_plus',
  'SIGIOBS(+)': '_refln.pdbx_I_plus_sigma',
  'IOBS(-)': '_refln.pdbx_I_minus',
  'SIGIOBS(-)': '_refln.pdbx_I_minus_sigma',
  # 'I-obs': '_refln.F_squared_meas',        # obsoleted
  # 'SIGI-obs': '_refln.F_squared_sigma',    # obsoleted
  'I-obs': '_refln.intensity_meas',
  'SIGI-obs': '_refln.intensity_sigma',
  'I-obs(+)': '_refln.pdbx_I_plus',
  'SIGI-obs(+)': '_refln.pdbx_I_plus_sigma',
  'I-obs(-)': '_refln.pdbx_I_minus',
  'SIGI-obs(-)': '_refln.pdbx_I_minus_sigma',
}

ccp4_to_cif_intensities = {
  # 'I': '_refln.F_squared_meas',            # which I to prefer?
  'I': '_refln.intensity_meas',              # This is used in PDB
  'SIGI': '_refln.intensity_sigma',          # This is used in PDB
  # 'SIGI': '_refln.F_squared_sigma',        # which SIGI to prefer?
  'I(+)': '_refln.pdbx_I_plus',
  'SIGI(+)': '_refln.pdbx_I_plus_sigma',
  'I(-)': '_refln.pdbx_I_minus',
  'SIGI(-)': '_refln.pdbx_I_minus_sigma',
  'DP': '_refln.pdbx_anom_difference',
  'SIGDP': '_refln.pdbx_anom_difference_sigma',
}

# -----------------------------------------------------------------------------
# amplitudes
phenix_to_cif_amplitudes = {
  'FOBS': '_refln.F_meas_au',
  'SIGFOBS': '_refln.F_meas_sigma_au',
  'FOBS(+)': '_refln.pdbx_F_plus',
  'SIGFOBS(+)': '_refln.pdbx_F_plus_sigma',
  'FOBS(-)': '_refln.pdbx_F_minus',
  'SIGFOBS(-)': '_refln.pdbx_F_minus_sigma',
  'F-obs': '_refln.F_meas_au',
  'SIGF-obs': '_refln.F_meas_sigma_au',
  'F-obs(+)': '_refln.pdbx_F_plus',
  'SIGF-obs(+)': '_refln.pdbx_F_plus_sigma',
  'F-obs(-)': '_refln.pdbx_F_minus',
  'SIGF-obs(-)': '_refln.pdbx_F_minus_sigma',
}

ccp4_to_cif_amplitudes = {
  'F': '_refln.F_meas_au',
  'SIGF': '_refln.F_meas_sigma_au',
  'FP': '_refln.F_meas_au',
  'SIGFP': '_refln.F_meas_sigma_au',
  'FC': '_refln.F_calc_au',
  'PHIC': '_refln.phase_calc',
  'PHIB': '_refln.phase_meas',
  'FOM': '_refln.fom',
  'FPART': '_refln.F_part_au',
  'PHIP': '_refln.phase_part',
  'F(+)': '_refln.pdbx_F_plus',
  'SIGF(+)': '_refln.pdbx_F_plus_sigma',
  'F(-)': '_refln.pdbx_F_minus',
  'SIGF(-)': '_refln.pdbx_F_minus_sigma',
  'DP': '_refln.pdbx_anom_difference',
  'SIGDP': '_refln.pdbx_anom_difference_sigma',
}

# -----------------------------------------------------------------------------
# map coefficients
phenix_to_cif_map_coefficients = {
  '2FOFCWT': '_refln.pdbx_FWT',
  'PH2FOFCWT': '_refln.pdbx_PHWT',
  'FOFCWT': '_refln.pdbx_DELFWT',
  'PHFOFCWT': '_refln.pdbx_DELPHWT',
}

ccp4_to_cif_map_coefficients = {
  'FWT': '_refln.pdbx_FWT',
  'PHWT': '_refln.pdbx_PHWT',
  'DELFWT': '_refln.pdbx_DELFWT',
  'DELPHWT': '_refln.pdbx_DELPHWT',
}

# -----------------------------------------------------------------------------
# Hendrickson-Lattman coefficients
phenix_to_cif_HL = {
  'HLA': '_refln.pdbx_HL_A_iso',
  'HLB': '_refln.pdbx_HL_B_iso',
  'HLC': '_refln.pdbx_HL_C_iso',
  'HLD': '_refln.pdbx_HL_D_iso',
}

ccp4_to_cif_HL = {
  'HLA': '_refln.pdbx_HL_A_iso',
  'HLB': '_refln.pdbx_HL_B_iso',
  'HLC': '_refln.pdbx_HL_C_iso',
  'HLD': '_refln.pdbx_HL_D_iso',
}

# -----------------------------------------------------------------------------
# R-free flag
phenix_to_cif_rfree = {
  'R-free-flags': '_refln.pdbx_r_free_flag'
}

ccp4_to_cif_rfree = {
  'FREE': '_refln.pdbx_r_free_flag',
}

# -----------------------------------------------------------------------------
# all mapping dictionaries
phenix_to_cif_labels_dict = dict()
phenix_to_cif_labels_dict.update(phenix_to_cif_intensities)
phenix_to_cif_labels_dict.update(phenix_to_cif_amplitudes)
phenix_to_cif_labels_dict.update(phenix_to_cif_map_coefficients)
phenix_to_cif_labels_dict.update(phenix_to_cif_HL)
phenix_to_cif_labels_dict.update(phenix_to_cif_rfree)

ccp4_to_cif_labels_dict = dict()
ccp4_to_cif_labels_dict.update(ccp4_to_cif_intensities)
ccp4_to_cif_labels_dict.update(ccp4_to_cif_amplitudes)
ccp4_to_cif_labels_dict.update(ccp4_to_cif_map_coefficients)
ccp4_to_cif_labels_dict.update(ccp4_to_cif_HL)
ccp4_to_cif_labels_dict.update(ccp4_to_cif_rfree)

# =============================================================================
# separate MTZ and CIF labels
mtz_intensity_labels = set(phenix_to_cif_intensities.keys())
mtz_intensity_labels.update(ccp4_to_cif_intensities.keys())

cif_intensity_labels = set(phenix_to_cif_intensities.values())
cif_intensity_labels.update(ccp4_to_cif_intensities.values())

mtz_amplitude_labels = set(phenix_to_cif_amplitudes.keys())
mtz_amplitude_labels.update(ccp4_to_cif_amplitudes.keys())

cif_amplitude_labels = set(phenix_to_cif_amplitudes.values())
cif_amplitude_labels.update(ccp4_to_cif_amplitudes.values())

mtz_map_coefficient_labels = set(phenix_to_cif_map_coefficients.keys())
mtz_map_coefficient_labels.update(ccp4_to_cif_map_coefficients.keys())

cif_map_coefficient_labels = set(phenix_to_cif_map_coefficients.values())
cif_map_coefficient_labels.update(ccp4_to_cif_map_coefficients.values())

mtz_HL_labels = set(phenix_to_cif_HL.keys())
mtz_HL_labels.update(ccp4_to_cif_HL.keys())

cif_HL_labels = set(phenix_to_cif_HL.values())
cif_HL_labels.update(ccp4_to_cif_HL.values())

mtz_rfree_labels = set(phenix_to_cif_rfree.keys())
mtz_rfree_labels.update(ccp4_to_cif_rfree.keys())

cif_rfree_labels = set(phenix_to_cif_rfree.values())
cif_rfree_labels.update(ccp4_to_cif_rfree.values())

# =============================================================================
# check new dictionary matches old dictionary (needs old version of
# mtz_as_cif.py)
# if (__name__ == '__main__'):
#   from iotbx.command_line import mtz_as_cif

#   assert (phenix_to_cif_labels_dict == mtz_as_cif.phenix_to_cif_labels_dict)
#   assert (ccp4_to_cif_labels_dict == mtz_as_cif.ccp4_to_cif_labels_dict)
def ccp4_label_from_cif(ciflabel):
  for mtzlabl,ciflabl in ccp4_to_cif_labels_dict.items():
    if ciflabel==ciflabl:
      return mtzlabl
  return None
