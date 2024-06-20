from ..cif.cif_browser import CifBrowserController
from ..filter import ComponentFilterObj, CompositeFilter, AtomFilterObj

class RestraintBrowserController(CifBrowserController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    
    #self.state.signals.ciffile_change.connect(self.update_file) # parent already implemented this
    self.state.signals.restraint_activated.connect(self.update_file) # should implement this

    # For restraints, the 'files' dropdown menu will only have one
    self.single_file_ref = None

  def update_file(self,ref):
    # Special callback to set the single file
    self.single_file_ref = ref

    super().update_file(ref)


  @property
  def cif_refs(self):
    return [self.single_file_ref]


  
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
        block_key = self.view.combobox_block.currentText()
        block_suffix = block_key.replace("_chem_comp_","")
        block_suffix_translator = {
          "bond":"bond",
          "angle":"angle",
          "tor":"dihedral",
          "chir":"chiral",
          "plane_atom":"plane"
        }
        block_suffix = block_suffix_translator[block_suffix]
        filter_obj_comp = ComponentFilterObj(name=comp_id,geometry_type=block_suffix)
        atom_columns = [col for col in df_sel.columns if "atom_id" in col]
        atom_filters = []
        for atom_col in atom_columns:

          atom_id = df_sel[atom_col].iloc[0]
          filter_obj_atom = AtomFilterObj(name=atom_id,geometry_type=block_suffix,mmcif_prefix=atom_col)
          atom_filters.append(filter_obj_atom)
        filter_obj = CompositeFilter([filter_obj_comp]+atom_filters)
        self.state.signals.geometry_filter_from_restraint.emit(filter_obj)
        print(dir(filter_obj))

        
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