from .atom_sites import AtomSites
from iotbx.data_manager import DataManager

class MolDF:
  """
  A molecule is composed of pandas dataframes
  """
  @classmethod
  def from_file_cif(cls,filename,**kwargs):
    atom_sites = AtomSites.from_file_cif(filename)
    return cls(atom_sites,params={"filename":filename},**kwargs)

  @classmethod
  def from_file_cif_via_iotbx(cls,filename,**kwargs):
    dm = DataManager()
    _ = dm.process_model_file(filename)
    model = dm.get_model()
    return cls.from_mmtbx_model(model,params={"filename":filename},**kwargs)
    
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

  @property
  def mmtbx_model(self):
    assert "mmtbx_model" in self.params, "No mmtbx model stored"
    return self.params["mmtbx_model"]

  # def select(self,arg,string_format='phenix',return_type="sites"):
  #   """
  #   An alias for self.atom_sites.select()
  #   """
  #   assert return_type in ["sites","query"], "Specify one of allowed return types"
  #   return self.atom_sites.select(arg,return_type=return_type)

    


    