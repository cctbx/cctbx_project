"""
The fundamental object to interact with the data/model/state of the GUI
"""
from __future__ import annotations # backwards compat string literal types
from pathlib import Path
import uuid
import hashlib
import json
from dataclasses import replace
from typing import Optional


from PySide2.QtCore import QObject, Signal
from mmtbx.model import manager as ModelManager


class Ref:
  """
    A Ref is the fundamental container for data. It provides:
      1. Decouple raw data from GUI-specific attributes
      2. Compositions of basic data structures that make sense in the GUI
      3. 'Instances' of data, where identical data appears in multiple places

  """
  _class_label_name = ""      # Subclass and provide a label
  EntryViewClass = None       # Associate an Entry View class (Entry is a widget class to display a Ref)
  EntryControllerClass = None # Associate an Entry Controller class (controller to control the entry view)

  def __init__(self,data: DataClassBase,show: bool = False):
    self._data = data
    self._uuid = self._generate_uuid()
    self._identifiers = {"uuid":self.uuid}
    self._label = None
    self._show = show
    self._entry = None # set later the entry controller object
    self.results = {} # program_name: result_ref
    self._active = False
    self._query = None



  # The 'data' is the fundamental data structure the ref is tracking. The data should be completely unaware
  #   it is inside a GUI. Meaning all GUI specific information related to that data should be stored 
  #   here in the Ref object
  @property
  def data(self):
    return self._data

  # Only one reference of a given type can be active at a time in the state. It is the implied
  #   subject for functions who don't get directly called, and need to ask the state who to operate on.
  #  
  @property
  def active(self):
    return self._active

  @active.setter
  def active(self,value):

    # toggle entry as active
    if self.entry:
      self.entry.active = value

    self._active = value

  # The show flag dictates whether a ref is presented in the GUI or not. 
  #   The controllers will use this to determine whether to load a ref into a gui widget or not

  @property
  def show(self):
    return self._show

  @ show.setter
  def show(self,value):
    self._show = value

  # Labels and Ids

  # Sometimes external programs will give data its own label/id. Meaning that a Ref may have multiple
  #   identifiers. 

  @property
  def uuid(self): # A truly unique id to label a Ref instance
    return self._uuid

  @property
  def identifiers(self): # A dictionary to access alternative ids by name, ie: identifiers["uuid"]
    return self._identifiers


  # Label is a text based label that the controllers will ask for when presenting a ref in the GUI
  @property
  def label(self):
    if self._label is None:
      if hasattr(self.data,"label") and self.data.label is not None:
        self._label = self.data.label

      elif hasattr(self.data,"filename") and self.data.filename is not None:
        try:
          # if key is a path, get the stem
          self._label = self._truncate_string(Path(self.data.filename).name)
        except:
          self._label = self._truncate_string(self.data.filename)
      else:
        self._label = self._class_label_name

    return self._label

  @label.setter
  def label(self,value):
    self._label = value
  
  # Utility methods
  def _truncate_string(self,path, max_len=20):
    # Truncate a long string to a specific number of characters,
    #   but show some of both the beginning and end. Useful for
    #   example when displaying absolute paths
    if len(path) > max_len:
      return path[:max_len // 2] + "..." + path[-max_len // 2:]
    else:
      return path


  @staticmethod
  def _generate_uuid(length: int=24):
    # Generate a truly unique identifier for a Ref instance
    full_uuid = str(uuid.uuid4())

    # Hash the UUID
    hashed_uuid = hashlib.sha1(full_uuid.encode()).hexdigest()

    # Truncate to the desired length
    short_uuid = hashed_uuid[:length]
    return short_uuid


# Subclasses of Ref appear below:

class ModelRef(Ref):
  """
  A Ref subclass for a molecular model. Specifically, a mmtbx.model.manager instance.

  The ModelRef can compose a GeometryRef and RestraintsRef to associate specific geometry and/or 
    generic restraints with the model.
  """
  _class_label_name = "model"
  def __init__(self,
    data: ModelManager, 
    geometry: Optional[GeometryRef] = None, 
    restraints: Optional[RestraintsRef] = None, 
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
  def model(self):
    return self.data.model

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



class SelectionRef(Ref):
  """
  A Ref subclass to track molecular selections. Selections can be generic, but a SelectionRef enforces
    an association between a Selection object and a specific molecular model
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