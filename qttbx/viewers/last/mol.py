from cctbx import crystal
from libtbx.utils import null_out
from .cif_io2 import CifInput
from .atom_sites2 import AtomSites
class DataFrameMol:
  """
  A mol is a composition of Pandas dataframes for each component (sites,bonds,etc)
  """


  @classmethod
  def from_mmcif_file(cls,file):
    cif_input = CifInput(file)
    sites = AtomSites.from_cif_input(cif_input)
    return cls(sites,cif_input=cif_input)

  @classmethod
  def from_mmtbx_model_via_mmcif(cls,model):
    cif_input = CifInput.from_mmtbx_model_via_mmcif(model)
    sites = AtomSites.from_cif_input(cif_input)
    cmol = cls(sites,cif_input=cif_input)
    cmol.model = model
    return cmol


  def __init__(self,sites,cif_input=None):
    assert isinstance(sites,AtomSites)
    self.cif_input = cif_input
    self._crystal_symmetry = None
    self._hierarchy = None
    self._sites = sites
    self._residues = None
    self._bonds = None
    self._angles = None
    self._dihedrals = None
    self._chirals = None
    self._planes = None
    self._grm = None
    self._model = None # mmtbx.model.manager


  # Component dataframes

  @property
  def hierarchy(self):
    return self._hierarchy

  @property
  def sites(self):
    return self._sites

  @sites.setter
  def sites(self,value):
    if not isinstance(value,AtomSites):
      value = AtomSites(value)
    self._sites = value

  @property
  def bonds(self):
    return self._bonds

  @bonds.setter
  def bonds(self,value):
    self._bonds = value

  @property
  def angles(self):
    return self._angles

  @angles.setter
  def angles(self,value):
    self._angles = value

  @property
  def dihedrals(self):
    return self._dihedrals

  @dihedrals.setter
  def dihedrals(self,value):
    self._dihedrals = value

  @property
  def chirals(self):
    return self._chirals

  @chirals.setter
  def chirals(self,value):
    self._chirals = value

  @property
  def planes(self):
    return self._planes

  @planes.setter
  def planes(self,value):
    self._planes = value

  # Misc properties
  @property
  def hierarchy(self):
    return self._hierarchy

  @property
  def crystal_symmetry(self):
    if self._crystal_symmetry is None:
      self.add_crystal_symmetry_if_necessary()
    return self._crystal_symmetry

  @property
  def grm(self):
    # geometry restraints manager
    return self._grm

  def add_crystal_symmetry_if_necessary(self,box_cushion=3):
    xyz = self.sites.xyz
    a,b,c = xyz.max(axis=0) - xyz.min(axis=0) + (2 * box_cushion)
    crystal_symmetry=crystal.symmetry((a,b,c, 90,90,90),1)
    self._crystal_symmetry = crystal_symmetry

  def build_hierarchy_if_necessary(self):
    if self._hierarchy is None:
      self._hierarchy = self.sites._build_hierarchy(sort=True)


  def add_i_seqs_from_hierarchy(self):
    self.build_hierarchy_if_necessary()
    i_seq_mapper = {atom.serial:atom.i_seq for atom in self.hierarchy.atoms()}
    for df_name in ["sites","bonds","angles","dihedrals",'chirals',"planes","nonbonded"]:
      df = getattr(self,df_name)
      if df is not None:
        for id_col in df.columns:
          if id_col.startswith("id"):
            if id_col == "id":
              iseq_col = "i_seq"
            elif "id_" in id_col:
              iseq_col = id_col.replace("id_","i_seq_")
            # add iseq column
            df[iseq_col] = df[id_col].map(i_seq_mapper).astype("Int64")

  # Functions for conversion

  @property
  def model(self):
    return self._model
  @model.setter
  def model(self,value):
    self._model = value

  def as_MMTBX_Model(self):
    from mmtbx.model import manager as ModelManager

    model = ModelManager(None,
                        pdb_hierarchy=self.hierarchy,
                        crystal_symmetry=self.crystal_symmetry,
                        #component_mol=self,
                        log=null_out())
    return model

