from typing import Optional

from mmtbx.model import manager as ModelManager
from qttbx.viewers.gui.model.refs import Ref

class ModelRef(Ref):
  """
  A Ref subclass for a molecular model. Specifically, a mmtbx.model.manager instance.

  The ModelRef can compose a GeometryRef and RestraintsRef to associate specific geometry and/or 
    generic restraints with the model.
  """
  _class_label_name = "model"
  def __init__(self,
    data: ModelManager, 
    geometry: Optional[object] = None, 
    restraints: Optional[object] = None, 
    show=True):
    super().__init__(data=data,show=show)
    self._geometry_ref = geometry
    self._restraints_ref = restraints


  @property
  def EntryControllerClass(self):
    from ..controller.models import ModelEntryController
    return ModelEntryController

  @property
  def EntryViewClass(self):
      from ..view.models import ModelEntryView
      return ModelEntryView

  @property
  def filename(self):
    model = self.model
    model_input = model.get_model_input()
    source_info = model_input.source_info()
    filename = None
    for line in source_info.split("\n"):
      if "file" in line:
        filename = line.split()[1]
        break
    return filename

  @property
  def model(self):
    return self.data

  # The optionally composed GeometryRef and RestraintsRef
  @property
  def geometry(self):
    return self._geometry_ref

  @geometry.setter
  def geometry(self,value):
    assert isinstance(value,(GeometryRef,type(None)))
    self._geometry_ref = value


  @property
  def has_geometry(self):
    return self.geometry is not None


  @property
  def restraints(self):
    return self._restraints_ref

  @restraints.setter
  def restraints(self,value):
    assert isinstance(value,(RestraintsRef,type(None)))
    self._restraints_ref = value

  @property
  def has_restraints(self):
    return self.restraints is not None

  def adopt_state(self,state):
    # An opportunity to do things when added to a state
    if self.__class__ not in state.active_refs:
      state.active_refs[self.__class__] = self