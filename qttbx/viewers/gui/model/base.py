"""
Base data structures which are managed by Ref instances. Data should be unaware of any other
objects in the GUI, and preferably immutable.
"""

import dataclasses
from dataclasses import dataclass, asdict, fields, make_dataclass,is_dataclass
import json
import uuid
import hashlib
from typing import Any, List, Tuple, Type, Optional
from collections import namedtuple
import pandas as pd

@dataclass(frozen=True)
class DataClassBase:


  def _coerce_types(self):
    for field_info in fields(self):
      value = getattr(self, field_info.name)
      if value is None:
        continue

      if field_info.type in [int, Optional[int]]:
        object.__setattr__(self, field_info.name, self._coerce_to_int(value))
      elif field_info.type in [float,Optional[float]]:
        object.__setattr__(self, field_info.name, self._coerce_to_float(value))
      elif field_info.type ==[str, Optional[str]]:
        object.__setattr__(self, field_info.name, self._coerce_to_str(value))


  @classmethod
  def from_defaults(cls):
    result = {}
    for f in fields(cls):
        # Check if the field has a default value
        if f.default is not dataclasses.MISSING:
            result[f.name] = f.default
    return result

  def to_dict(self):
    return self._to_dict(self)

  def _to_dict(self,obj):
    result = {}
    for field in fields(obj):
        value = getattr(obj, field.name)
        if is_dataclass(value):
            # Recursively convert nested dataclass to dict
            result[field.name] = self._to_dict(value)
        else:
            result[field.name] = value
    return result

  @classmethod
  def from_dict(cls,d,coerce_types=False):
    field_names = {f.name for f in fields(cls)}
    filtered_d = {k: v for k, v in d.items() if k in field_names}
    obj =  cls(**filtered_d)
    if coerce_types:
      obj._coerce_types()
    return obj

  def to_json(self, indent=None):
    return json.dumps(self.to_dict(),indent=indent)

  @classmethod
  def from_json(cls, json_str: str):
    d = json.loads(json_str)
    return cls.from_dict(d)



  @classmethod
  def from_phil_extract(cls,params,coerce_types=False):
    return cls.from_dict(params.__dict__,coerce_types=coerce_types)

  @property
  def hash(self):
    return f"{self.__class__.__name__}_{hash(self.to_json())}"

  @staticmethod
  def _generate_uuid(length: int=24):
    # Generate a UUID
    full_uuid = str(uuid.uuid4())

    # Hash the UUID
    hashed_uuid = hashlib.sha1(full_uuid.encode()).hexdigest()

    # Truncate to the desired length
    short_uuid = hashed_uuid[:length]
    return short_uuid

  # Functions to coerce to the right basic type
  def _coerce_to_int(self, value):
    try:
      return int(value)
    except (TypeError, ValueError):
      return value

  def _coerce_to_float(self, value):
    try:
      return float(value)
    except (TypeError, ValueError):
      return value

  def _coerce_to_str(self, value):
    return str(value)
