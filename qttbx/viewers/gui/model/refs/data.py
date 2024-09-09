from typing import Optional

from iotbx.data_manager import DataManager
from qttbx.viewers.gui.model.refs import Ref, ModelRef

class DataManagerRef(Ref):
  """
  A Ref subclass for a cctbx DataManager
  """
  _class_label_name = "data_manager"
  def __init__(self,
    data: DataManager):
    super().__init__(data=data)

  @property
  def data_manager(self):
    return self.data # alias

  def get_all_refs(self):
    # Extend this with maps/ etc
    return [self] + self.get_model_refs()

  def get_model_refs(self):
    # Get all model refs for models in data manager
    refs = []
    if hasattr(self.data_manager,"get_model_names"):
      for name in self.data_manager.get_model_names():
        model = self.data_manager.get_model(filename=name)
        ref = ModelRef(data=model)
        refs.append(ref)
    return refs