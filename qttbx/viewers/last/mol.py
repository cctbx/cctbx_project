from .atom_sites2 import AtomSites
from iotbx.data_manager import DataManager

class MolDF:
  """
  A molecule is composed of pandas dataframes
  """

  @classmethod
  def from_mmcif_file(cls,filename,**kwargs):
    sites = AtomSites.from_mmcif_file(filename=filename,**kwargs)
    return cls.from_atom_sites(sites,params={"filename":filename},**kwargs)

  @classmethod
  def from_atom_sites(cls,atom_sites,params={},**kwargs):
    return cls(atom_sites,params=params,**kwargs)

  @classmethod
  def from_mmtbx_model(cls,model,params={},**kwargs):
    atom_sites = AtomSites.from_mmtbx_model(model)
    params["mmtbx_model"] = model
    return cls(atom_sites,params=params,**kwargs)

  def __init__(self,atom_sites,params={},**kwargs):
    self._atom_sites = atom_sites
    self.params = params
    self.select = self.sites.select

  @property
  def sites(self):
    # alias of atom_sites
    return self.atom_sites

  @property
  def atom_sites(self):
    return self._atom_sites
