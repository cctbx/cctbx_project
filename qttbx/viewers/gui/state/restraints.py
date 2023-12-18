from dataclasses import dataclass
from typing import List, Optional


from .base import DataClassBase
from ...last.cctbx_utils import get_restraint_dfs_from_model

@dataclass(frozen=False)
class Restraint(DataClassBase):
  labels: List[str]
  i_seqs: List[int]



@dataclass(frozen=False)
class BondRestraint(Restraint):
  distance_ideal: float
  distance_model: float
  sigma: float
  delta: float
  residual: float
  weight: float
  slack: float
  sym_op_j: object
  rt_mx: object

@dataclass(frozen=False)
class AngleRestraint(Restraint):
  angle_ideal: float
  angle_model: float
  delta: float
  residual: float
  sigma: float
  weight: float


@dataclass(frozen=False)
class Restraints(DataClassBase):
  filename: Optional[str]
  bond: Optional[List[BondRestraint]]
  angle: Optional[List[AngleRestraint]]
  manager: Optional[object]
  components = {
    'bond': BondRestraint,
    'angle': AngleRestraint,
  }

  @classmethod
  def from_model_ref(cls,model_ref):
    model = model_ref.model
    dfs = get_restraint_dfs_from_model(model)
    instances_dict = {}
    for key,restraint_class in cls.components.items():
      df = dfs[key]
      instances = [restraint_class(**row._asdict()) for row in df.itertuples(index=False)]
      instances_dict[key] = instances

    return cls(filename=None,**instances_dict,manager=model.get_restraints_manager())
