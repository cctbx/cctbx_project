
from dataclasses import dataclass, asdict
import json

@dataclass
class DataClassBase:
  def to_dict(self):
    d = asdict(self)
    # TODO: Remove this when SelectionQuery is dataclass
    for key,value in list(d.items()):
      if hasattr(value,"to_dict"):
        d[key] = value.to_dict()
    return d
  
  def to_json(self, indent=None):
    return json.dumps(self.to_dict(),indent=indent)

  @classmethod
  def from_json(cls, json_str: str):
    params = json.loads(json_str)
    return cls(**params)