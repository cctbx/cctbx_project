
import dataclasses
from dataclasses import dataclass, asdict
import json
import uuid
import hashlib
from typing import Any, List, Tuple
from dataclasses import asdict, dataclass, fields
from typing import List, Type
import pandas as pd

import pandas as pd

@dataclass(frozen=True)
class DataClassBase:

  @staticmethod
  def defaults(cls):
    result = {}
    for f in fields(cls):
        # Check if the field has a default value
        if f.default is not dataclasses.MISSING:
            result[f.name] = f.default
    return result

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
    #d2["class_name"] = self.__class__.__name__
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
    field_names = {f.name for f in fields(cls)}
    filtered_d = {k: v for k, v in d.items() if k in field_names}
    return cls(**filtered_d)

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


@dataclass(frozen=True)
class ObjectRow(DataClassBase):
  pass


class ObjectFrame(list):
    def __init__(self, rows: List[Any], row_class: Type[Any]):
      super().__init__(rows)
      self.row_class = row_class

    def __repr__(self):
      return super().__repr__()

    def __iter__(self):
      return super().__iter__()

    def remove(self, value):
      assert isinstance(value,self.row_class), "Must call remove with an instance of the ObjectFrame's row class"
      hashes = self.hashes
      assert value.hash in hashes, "This value's hash does not appear in ObjectFrame's hashes"

      idx = hashes.index(value.hash)
      obj = self[idx]

      # Call the super class's remove method to perform the actual removal
      super().remove(obj)


    def iterrows(self):
      for row in self:
        yield row

    def iterfields(self, row: Any) -> Tuple[str, Any]:
      for field in fields(row):
        yield (field.name, getattr(row, field.name))
    @property
    def rows(self):
      return self

    @property
    def hash(self):
      return f"{self.__class__.__name__}_{self.row_class.__name__}_{hash(frozenset([row.hash for row in self]))}"
    @property
    def hashes(self):
      return [row.hash for row in self.rows]

    @property
    def df(self) -> pd.DataFrame:
      return self.to_df()

    def to_df(self) -> pd.DataFrame:
      rows_as_dicts = [asdict(row) for row in self]
      df = pd.DataFrame(rows_as_dicts)
      return df

    @classmethod
    def from_rows(cls, rows: List[Any]):
      types = [type(row) for row in rows]
      assert len(set(types)) == 1, "Mixed types cannot be used to init ObjectFrame"
      return cls(rows, row_class=type(rows[0]))

    @classmethod
    def from_df(cls, df: pd.DataFrame, row_class: Type[Any]) -> 'ObjectFrame':
      row_fields = {field.name for field in fields(row_class)}
      df_columns = set(df.columns)

      if not row_fields.issubset(df_columns):
        missing_fields = row_fields - df_columns
        extra_columns = df_columns - row_fields
        raise ValueError(f"DataFrame columns do not match ObjectRow fields. "
                          f"Missing fields: {missing_fields}, Extra columns: {extra_columns}")

      rows = [row_class(**{field: getattr(row, field) for field in row_fields}) for row in df.itertuples(index=False)]
      return cls(rows, row_class)
