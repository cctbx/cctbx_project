from typing import Optional

from qttbx.viewers.gui.model.refs import Ref, ModelRef
from qttbx.viewers.selection import Selection


class SelectionRef(Ref):
  """
  A Ref subclass to track molecular selections. 
  """
  _class_label_name = 'selection'
  def __init__(self,data: Selection,model: ModelRef,  show: Optional[bool] = True):
    assert show is not None, "Be explicit about whether to show selection ref in list of selections or not"
    assert ModelRef is not None, "Selection Ref cannot exist without reference model"
    assert data is not None, "Provide a Selection as the data"
    super().__init__(data=data, show=show)

    self._model_ref = model

  # The 'data' for a SelectionRef is a Selection
  @property
  def selection(self):
    return self.data

  @property
  def model(self):
    return self.model_ref.model

 
  @property
  def model_ref(self):
    return self._model_ref

  @property
  def number_of_atoms(self):
    return self.selection.bool.count(True)