'''
Canonical maps between DataManager datatypes and iotbx.file_reader.any_file
file types, plus the extension -> datatype table used by fast detection.

These live in iotbx.file_io (a low layer that does not import iotbx.data_manager)
so iotbx.data_manager can import them without a cycle. iotbx.data_manager
re-exports any_file_type / data_manager_type for backward compatibility.
'''

# datatype -> any_file file type
any_file_type = {
  'map_coefficients': 'hkl',
  'miller_array':     'hkl',
  'model':            'pdb',
  'ncs_spec':         'ncs',
  'phil':             'phil',
  'real_map':         'ccp4_map',
  'restraint':        'cif',
  'sequence':         'seq',
}

# any_file file type -> datatype (reverse), with the miller_array parent winning
# for hkl and the legacy mtz token accepted
data_manager_type = {value: key for key, value in any_file_type.items()}
data_manager_type['hkl'] = 'miller_array'
data_manager_type['mtz'] = 'miller_array'

# Sentinel: a .cif/.mmcif extension is ambiguous (model/restraint/miller_array)
# and is resolved by an authoritative xcif parse of the content, not the
# extension alone (see detection._cif_datatypes).
CIF_SENTINEL = 'cif?'

# file extension (lowercase, no dot) -> DataManager datatype
extension_to_datatype = {
  'pdb': 'model', 'ent': 'model',
  'mtz': 'miller_array', 'hkl': 'miller_array', 'sca': 'miller_array',
  'cns': 'miller_array', 'cv': 'miller_array', 'ref': 'miller_array',
  'fobs': 'miller_array',
  'cif': CIF_SENTINEL, 'mmcif': CIF_SENTINEL,
  'ccp4': 'real_map', 'mrc': 'real_map', 'map': 'real_map',
  'ncs': 'ncs_spec', 'ncs_spec': 'ncs_spec',
  'params': 'phil', 'eff': 'phil', 'def': 'phil', 'phil': 'phil',
  'param': 'phil',
  'fa': 'sequence', 'faa': 'sequence', 'seq': 'sequence', 'pir': 'sequence',
  'dat': 'sequence', 'fasta': 'sequence',
  'json': 'json',
  'yaml': 'yaml', 'yml': 'yaml',
}
