"""
The fundamental object to interact with the data/model/state of the GUI
"""
from pathlib import Path
import uuid
import hashlib
import json
from typing import Optional


from PySide2.QtCore import QObject, Signal


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

  def __init__(self,data: object,show: bool = False):
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
  # @property
  # def active(self):
  #   return self._active

  # @active.setter
  # def active(self,value):

  #   # toggle entry as active
  #   if self.entry:
  #     self.entry.active = value

  #   self._active = value

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


  def adopt_state(self,state):
    # An opportunity to do things when added to a state
    pass
