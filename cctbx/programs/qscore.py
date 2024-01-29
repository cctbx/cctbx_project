from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
from libtbx import group_args
from iotbx.pdb import input as pdb_input
from mmtbx.model import manager as model_manager
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
from cctbx.maptbx import qscore
import numpy as np

# =============================================================================

class Program(ProgramTemplate):

    description = """
  Perform a Qscore analysis for map-model fit
  """

    datatypes = ['phil', 'model', 'real_map']

    master_phil_str = """
    include scope cctbx.maptbx.qscore.master_phil_str
    """

    def validate(self):
        pass

    def run(self):
        print("Running")

        # do selection
        mmm = self.data_manager.get_map_model_manager()
        if self.params.qscore.selection != None:
            selection = np.where(mmm.model().selection(self.params.qscore.selection).as_numpy_array())[0]
            if len(selection)==0:
                print("Finished... nothing selected")
                self.result = group_args()
                return
            self.params.qscore.selection = selection

        # calculate shells

        if len(self.params.qscore.shells) ==0 :
           start = self.params.qscore.shell_radius_start
           stop = self.params.qscore.shell_radius_stop
           num = self.params.qscore.shell_radius_num
           shells = list(np.linspace(
            start,
            stop,
            num,
            endpoint=True))
           for shell in shells:
              self.params.qscore.shells.append(shell)


        # unsort model, make new mmm
        model = self._unsort_model(mmm.model())
        # ignore hydrogens
        model = model.select(model.selection("not element H"))

        # make mmm
        mmm.set_model(model,overwrite=True)


        # run qscore
        qscore_result= qscore.calc_qscore(
            mmm,
            self.params.qscore,
            log=self.logger)


        self.result = group_args(**qscore_result)



    def get_results(self):
        return self.result


    def _unsort_model(self,model):

      # A very hacky way to get a model with atoms unsorted. Only works for pdb
      # TODO: Need to add sort kwarg to data_manager.process_model_file somehow
      model_input = model.get_model_input()
      if 'cif' not in str(model_input.__class__):
        h = pdb_input(source_info=None,lines=model_input.as_pdb_string()).construct_hierarchy(sort_atoms=False)
        atom_name_list = [s.strip() for s in h.atoms().extract_name()]
        scatterer_list = [s.strip() for s in h.atoms().extract_element()]

        model = model_manager.from_sites_cart(h.atoms().extract_xyz(),
                                                    crystal_symmetry=model.crystal_symmetry(),
                                                    atom_name_list = atom_name_list,
                                                    scatterer_list=scatterer_list
                                                    )
      return model

    @staticmethod
    def _mmm_from_model_and_map(model,map):
      map_data = map.map_data()
      mm = map_manager(map_data=map_data,
                unit_cell_grid=map_data.accessor().all(),
                unit_cell_crystal_symmetry=model.unit_cell_crystal_symmetry(),
                wrapping=False)
      mmm = map_model_manager(model=model,map_manager=mm)
      return mmm
