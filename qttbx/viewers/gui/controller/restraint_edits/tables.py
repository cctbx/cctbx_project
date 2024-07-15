

from ..table import TableController

from ...state.table import  PandasTableModel
from ...state.edits import (
  EditData,
  BondEdit,
  AngleEdit,
  DihedralEdit,
)



class EditsTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  row_class = EditData # Generic, use subclasses
  restraint_name =  None # Subclass and fill this in (bond,angle,etc)
  #columns_to_include= ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  # supress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old", "labels_compositional"]
  # rename_columns = {"ideal_new":"Ideal",
  #                   "sigma_new":"Sigma",
  #                   "action":"Action",
  #                   "atom_id_1":"Atom 1",
  #                   "atom_id_2":"Atom 2"}

  # column_prefixes_to_include = ["atom_id"]

  @staticmethod
  def transform_to_dict(nested_list,prefix="Label"):
    # Transform a list of lists to a dictionary of lists
    # Determine the number of sublists
      num_sublists = len(nested_list)
      # Determine the length of each sublist
      sublist_length = len(nested_list[0])
      
      # Initialize an empty dictionary to store the results
      result_dict = {}
      
      # Iterate over the indices of the sublists
      for i in range(sublist_length):
        # Create a label for each index
        label = f'{prefix}_{i+1}'
        # Gather all elements at index i from each sublist
        result_dict[label] = [nested_list[j][i] for j in range(num_sublists)]
      
      return result_dict


  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.edits_ref = None

    # Signals
    self.state.signals.edits_added.connect(self.update)
    self.view.table_view.removeEdit.connect(self.remove_edit)
    self.view.write_button.clicked.connect(self.write_edits)

  def update(self,edits_ref):
    if edits_ref.data.row_class==self.row_class:
      self.edits_ref = edits_ref
      self.dataframe = edits_ref.data.df
      if len(self.dataframe)>0:
        # unpack labels
        label_dict = self.transform_to_dict(self.dataframe["labels_compositional"].tolist())
        for key,value in label_dict.items():
          self.dataframe[key] = value
        
        self.table_model = PandasTableModel(self.dataframe,display_columns=self.display_columns)
        self.parent.view.toggle_tab_visible(self.title,show=True)
      else:
        self.state.signals.remove_ref.emit(edits_ref)
        self.parent.view.toggle_tab_visible(self.title,show=False)

  def remove_edit(self,row_dict):
    if self.edits_ref is not None:
      edit = self.row_class.from_dict(row_dict)
      self.edits_ref.data.remove(edit)
      self.update(self.edits_ref)

  def write_edits(self):
    self.parent.write_all_edits()


# TODO: Rename to include 'edit' name
class BondTableController(EditsTableTabController):
  title = "Bonds"
  restraint_name = "bond"
  row_class = BondEdit
  display_columns = ["action","ideal_old","ideal_new","sigma_old","sigma_new","Label"]

class AngleTableController(EditsTableTabController):
  title = "Angles"
  restraint_name = "angle"
  row_class = AngleEdit
  display_columns = ["action","ideal_old","ideal_new","sigma_old,sigma_new","Label"]


class DihedralTableController(EditsTableTabController):
  title = "Dihedrals"
  restraint_name = "dihedral"
  row_class = DihedralEdit

class ChiralTableController(EditsTableTabController):
  title = "Chirals"
  restraint_name = "chirality"

class PlanarityTableController(EditsTableTabController):
  title = "Planes"
  restraint_name = "plane"

class NonbondedTableController(EditsTableTabController):
  title = "Nonbonded"
  restraint_name = "nonbonded"



