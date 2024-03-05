
from dataclasses import dataclass, asdict
import json
import uuid
import hashlib
from typing import Optional

@dataclass
class DataClassBase:
  def update(self, **kwargs):
    for key, value in kwargs.items():
      if hasattr(self, key):
        if value is not None:
          setattr(self, key, value)


  def __post_init__(self):
    pass
  def to_dict(self):
    """"
    Conversion with to_dict (and so json) does not automatically
    include everything. If a full dict is desired use asdict directly.
    """
    d1 = asdict(self)
    d2 = {}
    d2["class_name"] = self.__class__.__name__
    # TODO: Remove this when SelectionQuery is dataclass
    for key,value in list(d1.items()):
      if hasattr(value,"to_dict"):
        d2[key] = value.to_dict()
      elif isinstance(value,(str,int,float,dict,list,tuple)):
        d2[key] = value
      else:
        pass # value may not be json serializable.
    return d2

  def to_json(self, indent=None):
    return json.dumps(self.to_dict(),indent=indent)

  @classmethod
  def from_json(cls, json_str: str):
    d = json.loads(json_str)
    return cls.from_dict(d)

  @classmethod
  def from_dict(cls,d):
    if 'class_name' in d:
      del d["class_name"]
    return cls(**d)

  @staticmethod
  def _generate_uuid(length: int=24):
    # Generate a UUID
    full_uuid = str(uuid.uuid4())

    # Hash the UUID
    hashed_uuid = hashlib.sha1(full_uuid.encode()).hexdigest()

    # Truncate to the desired length
    short_uuid = hashed_uuid[:length]
    return short_uuid