"""
Data structures to descrine the internal state of the molstar viewer, but in python.
"""
from dataclasses import dataclass
from typing import List, Optional, Dict

from .base import DataClassBase


@dataclass(frozen=True)
class Component(DataClassBase):
  phenixKey: str
  key: str
  representations: List[str]


@dataclass(frozen=True)
class Representation(DataClassBase):
  phenixKey: str
  name: str


@dataclass(frozen=True)
class Structure(DataClassBase):
  phenixReferenceKey: str
  phenixKey: str
  data_id: str
  key: str
  components: List[Component]



@dataclass(frozen=True)
class Reference(DataClassBase):
  id_viewer: Optional[str]
  id_molstar: Optional[str]
  structures: List[Structure]

  @classmethod
  def from_default(cls):
    return cls(
      id_viewer = None,
      id_molstar = None,
      structure = None
    )

  @property
  def representations(self):
    output = []
    for structure in self.structures:
      for component in structure.components:
        for repr_name in component.representations:
          output.append(repr_name)
    return output


@dataclass(frozen=True)
class MolstarState(DataClassBase):
  has_synced: bool
  references: Dict[str,Reference]

  @classmethod
  def from_empty(cls):
    return cls(references={},has_synced=False)

  @classmethod
  def from_dict(cls,state_dict):
    has_synced = state_dict["has_synced"]
    phenix_state =  cls(references = {}, has_synced=has_synced)
    # update state from dict
    for ref_dict in state_dict["references"]:
      id_molstar = ref_dict["molstarKey"]
      id_viewer = ref_dict["phenixKey"]

      ref = Reference(id_molstar=id_molstar, id_viewer = id_viewer, structures = [])
      phenix_state.references[id_viewer] = ref
      # now modify within a ref
      for structure_dict in ref_dict['structures']:
        structure = Structure(
          phenixKey=structure_dict["phenixKey"],
          phenixReferenceKey=structure_dict['phenixReferenceKey'],
          data_id=structure_dict["data_id"],
          key=structure_dict["key"],
          components=[])
        ref.structures.append(structure)
        for component_dict in structure_dict['components']:
          component = Component(phenixKey=component_dict["phenixKey"],representations=[],key=component_dict["key"])
          structure.components.append(component)

          for representation_dict in component_dict["representations"]:
            representation = Representation(phenixKey=representation_dict["phenixKey"],name=representation_dict["name"])

            component.representations.append(representation)

    return phenix_state