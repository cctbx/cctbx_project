

import dataclasses
from dataclasses import dataclass, asdict, fields, make_dataclass
import json
import uuid
import hashlib
from typing import Any, List, Tuple, Type, Optional
from collections import namedtuple
import pandas as pd

@dataclass(frozen=True)
class DataClassBase:


  def update(self, **kwargs):
    for key, value in kwargs.items():
      if hasattr(self, key):
        if value is not None:
          setattr(self, key, value)

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
  def defaults(cls):
    result = {}
    for f in fields(cls):
        # Check if the field has a default value
        if f.default is not dataclasses.MISSING:
            result[f.name] = f.default
    return result

  def _as_dict(self):
    return asdict(self)

  def to_dict(self):
    """"
    Conversion with to_dict (and so json) does not automatically
    include everything. If a full dict is desired use asdict directly.
    """
    d1 = self._asdict(self)
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
  def from_dict(cls,d,coerce_types=False):
    field_names = {f.name for f in fields(cls)}
    filtered_d = {k: v for k, v in d.items() if k in field_names}
    obj =  cls(**filtered_d)
    if coerce_types:
      obj._coerce_types()
    return obj

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


@dataclass(frozen=True)
class ObjectRow(DataClassBase):
  pass


class ObjectFrame(list):
    def __init__(self, rows: List[Any] = [], row_class: Type[Any] = DataClassBase):
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

    def to_df(self, unstack: bool = True) -> pd.DataFrame:
      rows_as_dicts = [asdict(row) for row in self]
      df = pd.DataFrame(rows_as_dicts)
      if unstack:
        if "i_seqs" in df.columns:
          df = self.unpack_list_column(df, 'i_seqs',remove_plural_s=True)
      return df

    @staticmethod
    def unpack_list_column(df: pd.DataFrame, col_name: str,remove_plural_s=True) -> pd.DataFrame:
      max_len = df[col_name].apply(lambda x: len(x) if isinstance(x, list) else 0).max()
      for i in range(max_len):
        if remove_plural_s and col_name.endswith("s"):
          col_name_new = col_name[:-1]
        else:
          col_name_new = col_name
        df[f'{col_name_new}_{i+1}'] = df[col_name].apply(lambda x: x[i] if isinstance(x, list) and i < len(x) else None)
      return df

    @staticmethod
    def collapse_cols(df,prefix):
      cols = [col for col in df.columns if col.startswith(prefix)]
      new_col = prefix+"s"
      df[new_col] = df.apply(lambda row: [row[col] for col in cols], axis=1)
      return df

    @classmethod
    def from_rows(cls, rows: List[Any]):
      types = [type(row) for row in rows]
      assert len(set(types)) == 1, "Mixed types cannot be used to init ObjectFrame"
      return cls(rows, row_class=type(rows[0]))

    @classmethod
    def from_df(cls, df: pd.DataFrame, row_class: Optional[Type[DataClassBase]] = None, name: str = "objects") -> 'ObjectFrame':
      if row_class:
        row_fields = {field.name for field in fields(row_class)}
        df_columns = set(df.columns)

        # Check for columns that start with 'i_seq_' and combine them into a single 'i_seqs' column
        i_seq_columns = [col for col in df_columns if col.startswith('i_seq_')]
        if i_seq_columns:
          df['i_seqs'] = df[i_seq_columns].apply(lambda row: [val for val in row if pd.notna(val)], axis=1)
          df.drop(columns=i_seq_columns, inplace=True)

        # Ensure DataFrame has all required columns by adding missing ones with default values
        for field in row_fields:
          if field not in df.columns:
            df[field] = None  # or some default value appropriate for your use case

        # Filter the DataFrame to only include the necessary columns
        df = df[[col for col in df.columns if col in row_fields]]

        if not row_fields.issubset(set(df.columns)):
          missing_fields = row_fields - set(df.columns)
          raise ValueError(f"DataFrame columns do not match ObjectRow fields. Missing fields: {missing_fields}")

        rows = []
        for row in df.itertuples(index=False, name='Pandas'):
          row_dict = {field: getattr(row, field) for field in row_fields}
          rows.append(row_class(**row_dict))
      else:
        # Check for columns that start with 'i_seq_' and combine them into a single 'i_seqs' column
        df_columns = set(df.columns)
        i_seq_columns = [col for col in df_columns if col.startswith('i_seq_')]
        if i_seq_columns:
          df['i_seqs'] = df[i_seq_columns].apply(lambda row: [val for val in row if pd.notna(val)], axis=1)
          df.drop(columns=i_seq_columns, inplace=True)

        # Dynamically generate fields
        RecordNT = namedtuple(name, df.columns)
        named_tuples = [RecordNT(*row) for row in df.itertuples(index=False, name=None)]

        # Dynamically create a dataclass based on the fields in the named tuple and subclass DataClassBase
        fields_list = [(col, List[int] if col == 'i_seqs' else type(val)) for col, val in zip(df.columns, df.iloc[0])]
        row_class = make_dataclass('GenericRow', fields_list, bases=(DataClassBase,),frozen=True)

        # Convert named tuples to dataclass instances
        rows = [row_class(*nt) for nt in named_tuples]

      return cls(rows, row_class)
