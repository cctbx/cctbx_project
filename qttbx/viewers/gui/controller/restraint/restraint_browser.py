from pathlib import Path

from PySide2.QtWidgets import QFileDialog

from ...state.ref import RestraintsRef
from ..cif.cif_browser import CifBrowserController
from ..filter import ComponentFilterObj, CompositeFilter, AtomFilterObj

class RestraintBrowserController(CifBrowserController):
  """
  Inherits from TableController -> CifBrowserController -> RestraintBrowserController
  """
  # Modifications to data blocks
  rename_data_blocks = {}
  rename_data_items = {"_chem_comp_bond":"Bonds",
                      "_chem_comp_angle":"Angles",
                      "_chem_comp_tor":"Dihedrals",
                      "_chem_comp_chir":"Chirals",
                      "_chem_comp_plane_atom":"Planes"}
  suppress_data_blocks = ["data_comp_list"]
  suppress_data_items = ["_chem_comp_tree","_chem_comp_atom","_chem_comp_rotamer_info"]

  # Modifications to columns
  display_columns = [] 
  column_display_names = {
    "value_dist":"Ideal",
    "value_dist_neutron":"Ideal neutron",
    "value_angle":"Ideal",
    "value_dist_esd":"Sigma",
    "value_angle_esd":"Sigma",
    "dist_esd":"Sigma"
  
  }
  editable_columns = ["value_dist","value_dist_neutron","value_dist_esd", #bonds
                      "value_angle", "value_angle_esd", # angles
                       "period",# Dihedrals
                       "volume_sign", # Chirals
                       "dist_esd", # planes
                      ] # if not empty, enforced as exclusive
  non_editable_columns = [] # excluded regardless of presence in editable list
  
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    
    #self.state.signals.ciffile_change.connect(self.update_file) # parent already implemented this
    #self.state.signals.restraint_activated.connect(self.update_file) # should implement this
    self.state.signals.restraints_change.connect(self.update_file)



  def update_file(self,ref):
    # incoming ref will be RestraintsRef not RestraintRef. Pick first one
    if isinstance(ref,RestraintsRef):
      ref = ref.restraints[0]

    super().update_file(ref)


  @property
  def cif_refs(self):
    if self.state.active_model_ref:
      return self.state.active_model_ref.restraints_ref.restraints
    else:
      return []

  @property
  def filenames(self):
    # implemented to enable modification in subclasses
    filenames = [ref.data.filename.replace("data_","").replace(".cif","") for ref in self.cif_refs]
    return filenames

  def save(self,*args):
    # Opens a save file dialog and returns the selected file path and filter

    path = Path(self.cif_ref.data.filepath)
    suggested_path = str(Path(Path.home(),path.stem+"_edited"+"".join(path.suffixes)))
    filepath, _ = QFileDialog.getSaveFileName(self.view, "Save File", suggested_path, "All Files (*);")

    self.rename_df_dict_undo()
    if filepath:
        print(f"File selected for saving: {filepath}")
        write_dataframes_to_cif_file(self.df_dict,str(Path(filepath)))
  
  def on_selection_changed(self, selected, deselected):
    # Send a selection signal out to focus in graphics
    df_sel = self.view.table_view.selected_rows()

    # Override clicking on row
    # Apply a filter over in geometry tab
    if len(df_sel)==1:
      if "comp_id" in df_sel.columns and "atom_id_1" in df_sel.columns:
        # a candidate for filtering
        print("Will try to filter by restraint")
        comp_id = df_sel["comp_id"].iloc[0]
        item_key = self.view.combobox_data_item.currentText()
        item_suffix = item_key.replace("_chem_comp_","")
        item_suffix_translator = {
          "bond":"bond",
          "Bonds":"bond",
          "angle":"angle",
          "Angles":"angle",
          "tor":"dihedral",
          "Dihedrals":"dihedral",
          "chir":"chiral",
          "Chirals":"chiral",
          "plane_atom":"plane",
          "Planes":"plane"
        }
        item_suffix = item_suffix_translator[item_suffix]
        filter_obj_comp = ComponentFilterObj(name=comp_id,geometry_type=item_suffix)
        atom_columns = [col for col in df_sel.columns if "atom_id" in col]
        atom_filters = []
        for atom_col in atom_columns:

          atom_id = df_sel[atom_col].iloc[0]
          filter_obj_atom = AtomFilterObj(name=atom_id,geometry_type=item_suffix,mmcif_prefix=atom_col)
          atom_filters.append(filter_obj_atom)
        filter_obj = CompositeFilter([filter_obj_comp]+atom_filters)
        self.state.signals.geometry_filter_from_restraint.emit(filter_obj)
        
      else:
        print("Not a candidate row for filtering by restraint")



    # if len(df_sel)==1:
    #   if df_sel is not None:
    #     # switch to atom picking level
    #     self.state.signals.picking_level.emit("atom")
    #     flattened_i_seqs = [item for sublist in df_sel['i_seqs'] for item in sublist]
    #     flattened_i_seqs = [int(e) for e in flattened_i_seqs if pd.notna(e)]
    #     selection = Selection.from_i_seqs(self.state.mol.sites,flattened_i_seqs)
    #     ref = SelectionRef(data=selection,model_ref=self.state.active_model_ref,show=False)
    #     self.state.add_ref(ref)
    #     self.state.active_selection_ref = ref
    #   else:
    #     print("no atoms returned as query")