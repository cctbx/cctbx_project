from iotbx import pdb
import iotbx.pdb.parser
import iotbx.pdb.cryst1_interpretation
import iotbx.pdb.atom

class empty: pass

def clean_model_serial_list(model_serial_list):
  result = []
  for model_serial in model_serial_list:
    if (model_serial in result):
      model_serial = max(result) + 1
      if (model_serial in model_serial_list):
        model_serial = max(model_serial_list) + 1
    result.append(model_serial)
  return result

class stage_1:

  def __init__(self, file_name=None, raw_records=None):
    assert [file_name, raw_records].count(None) == 1
    if (raw_records is None):
      raw_records = open(file_name)
    columns_73_76_eval = pdb.parser.columns_73_76_evaluator(
      raw_records=raw_records)
    raw_records = columns_73_76_eval.raw_records
    self.ignore_columns_73_and_following = columns_73_76_eval.is_old_style
    self.crystal_symmetry = None
    self.remark_290_records = []
    self.ter_indices = []
    self.break_indices = []
    self.atom_attributes_list = []
    self.conect_records = []
    self.link_records = []
    self.ssbond_records = []
    self.sltbrg_records = []
    model_serial = None
    model_serial_list = []
    state = empty()
    self.state = state
    for state.line_number,state.raw_record in enumerate(raw_records):
      record_name = state.raw_record[:6]
      if (record_name == "CRYST1"):
        self.crystal_symmetry = pdb.cryst1_interpretation.crystal_symmetry(
          cryst1_record=state.raw_record,
          line_number=state.line_number)
      elif (state.raw_record.startswith("REMARK 290 ")):
        self.remark_290_records.append(state.raw_record.rstrip())
      elif (record_name == "MODEL "):
        model_serial = self.parse_record().serial
        model_serial_list.append(model_serial)
      elif (record_name == "ENDMDL"):
        model_serial = None
      elif (record_name.rstrip() == "TER"):
        self.ter_indices.append(len(self.atom_attributes_list))
      elif (record_name.rstrip() == "BREAK"):
        self.break_indices.append(len(self.atom_attributes_list))
      elif (record_name in ("ATOM  ", "HETATM")):
        atom_attributes = pdb.atom.attributes().set_from_ATOM_record(
          self.parse_record())
        if (model_serial is None):
          atom_attributes.MODELserial = -1
        else:
          atom_attributes.MODELserial = len(model_serial_list)-1
        self.atom_attributes_list.append(atom_attributes)
      elif (record_name == "SIGATM"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_SIGATM_record(
            self.parse_record())
      elif (record_name == "ANISOU"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_ANISOU_record(
            self.parse_record())
      elif (record_name == "SIGUIJ"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_SIGUIJ_record(
            self.parse_record())
      elif (record_name == "CONECT"):
        self.conect_records.append(state.raw_record)
      elif (record_name == "LINK  "):
        self.link_records.append(state.raw_record)
      elif (record_name == "SSBOND"):
        self.ssbond_records.append(state.raw_record)
      elif (record_name == "SLTBRG"):
        self.sltbrg_records.append(state.raw_record)
    del self.state
    self.model_serial_list = clean_model_serial_list(model_serial_list)
    if (len(self.model_serial_list) > 0):
      for atom_attributes in self.atom_attributes_list:
        if (atom_attributes.MODELserial < 0): continue
        atom_attributes.MODELserial \
          = self.model_serial_list[atom_attributes.MODELserial]
    self._selection_cache = None

  def parse_record(self):
    return pdb.parser.pdb_record(
      raw_record=self.state.raw_record,
      line_number=self.state.line_number,
      ignore_columns_73_and_following=self.ignore_columns_73_and_following)

  def selection_cache(self):
    if (self._selection_cache is None):
      self._selection_cache = pdb.atom.selection_cache(
        atom_attributes_list=self.atom_attributes_list)
    return self._selection_cache
