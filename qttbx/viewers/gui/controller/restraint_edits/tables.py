

from ..table import TableController

from ...state.table import  PandasTableModel
from ...state.edits import (
  EditData,
  BondEdit,
  AngleEdit,
  DihedralEdit,
  ChiralEdit,
  PlaneEdit,
)



class EditsTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  row_class = EditData # Generic, use subclasses
  restraint_name =  None # Subclass and fill this in (bond,angle,etc)
  display_columns = []
  rename_columns = {}

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
    self.view.read_button.clicked.connect(self.read_edits)
    self.view.write_button.clicked.connect(self.write_edits)
    

  def update(self,edits_ref):
    if edits_ref.data.row_class==self.row_class:
      self.edits_ref = edits_ref
      self.dataframe = edits_ref.data.df
      if len(self.dataframe)>0:        
        self.table_model = PandasTableModel(self.dataframe,
        display_columns=self.display_columns,
        rename_columns = self.rename_columns,
        capitalize=True,
        remove_underscores = True,
        transpose=False)
        self.parent.view.toggle_tab_visible(self.title,show=True)

  def remove_edit(self,row_dict):
    if self.edits_ref is not None:
      edit = self.row_class.from_dict(row_dict)
      self.edits_ref.data.remove(edit)
      self.update(self.edits_ref)

  def read_edits(self):
    # This is not good practice to reach over to another tab, but works
    self.parent.parent.molstar.viewer_controls.open_edits_file_dialog()

  def write_edits(self):
    self.parent.write_all_edits()

  def on_selection_changed(self,*args):
    pass

# TODO: Rename to include 'edit' name
class BondTableController(EditsTableTabController):
  title = "Bonds"
  restraint_name = "bond"
  row_class = BondEdit
  display_columns = ["action","atom_label","distance_ideal","sigma",]

class AngleTableController(EditsTableTabController):
  title = "Angles"
  restraint_name = "angle"
  row_class = AngleEdit
  display_columns = ["action","atom_label","angle_ideal","sigma",]


class DihedralTableController(EditsTableTabController):
  title = "Dihedrals"
  restraint_name = "dihedral"
  row_class = DihedralEdit
  display_columns = ["action","atom_label","angle_ideal","sigma","harmonic"]


class ChiralTableController(EditsTableTabController):
  title = "Chirals"
  restraint_name = "chirality"
  row_class = ChiralEdit
  display_columns = ["action","ideal_old","ideal_new","sigma_old","sigma_new","Label"]


class PlanarityTableController(EditsTableTabController):
  title = "Planes"
  restraint_name = "plane"

class NonbondedTableController(EditsTableTabController):
  title = "Nonbonded"
  restraint_name = "nonbonded"



