import random
from io import StringIO
from libtbx import group_args
from itertools import chain
import sys
import time

from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QTimer

from qttbx.viewers.gui.model.selection import Selection, PhenixParser
from qttbx.viewers.gui.model.state import State
from qttbx.viewers.gui.controller.molstar import MolstarController
from qttbx.viewers.gui.view.molstar import *

class SelectionTester:
    
    test_params_default = {
        "operators":["=",">", "<"],
        "bfactors":[20.0, 40.0, 80.0, 100.0],
        "occupancies":[.20, .40 .80, 1.0],
        "logic_tests":["and"] + ["or"]*4 + ["and not"], # The logic tests to draw from
        "expression_sizes":[1]*5 + [2]*5 + [3]*3 + [4]*2 + [5,6], # count of expressions to combine for a test,
        "scale_n_tests_by": 100, # scale factor for generating a quantity of random tests
    }

    @classmethod
    def from_file(cls,filename):
        dm = DataManager()
        dm.process_model_file(str(filename))
        test_params = cls.test_params_default
        return cls(dm.get_model(),test_params)
        
    def __init__(self,model,test_params=None):
        assert model
        self.model = model
        if not test_params:
            test_params = self.test_params_default
            
        chains, resnames, resseqs, names, alts = self.generate_options(self.model)
        self.chains = chains
        self.resnames = resnames
        self.resseqs = resseqs
        self.names = names
        self.alts = alts
        out = StringIO()
        self.model.composition().show(out)
        self.log(out.getvalue())
        self.log("Altlocs: ",self.alts)
        for key,value in test_params.items():
            setattr(self,key,value)

        self.selection_fragments = None
        self.valid_fragments = None
        self.selection_expressions = None
        self.valid_expressions = None
        self.failures = 0

    def log(self,*args):
        print(*args)
    
    def generate_options(self,model):
        chain_options = sorted(list(set(model.chain_ids())))
        resname_options = sorted({
          resname.strip()
          for chain in model.chains()
          for residue_group in chain.residue_groups()
          for resname in residue_group.unique_resnames()
        })
        resseq_options = sorted({
          residue_group.resseq_as_int()
          for chain in model.chains()
          for residue_group in chain.residue_groups()
        })
        resseq_options = sorted([int(e) for e in resseq_options ])
        resseq_options = [str(e) for e in resseq_options]
        atom_options = sorted(list(set([a.name.strip() for a in model.get_atoms()])))
        alts_options = sorted(list(set([a.parent().altloc for a in model.get_atoms()])))
        return chain_options, resname_options, resseq_options, atom_options, alts_options
    
    def operator_choose(self):
        return random.choice(self.operators)
        
    def chain_choose(self):
        return random.choice(self.chains)
    
    def resname_choose(self):
        return random.choice(self.resnames)
        
    def resseq_choose(self):
      space = list(range(0,len(self.resseqs)))
      low = random.choice(space)
      high = random.choice(space[low:])
      return low, high
    
    def names_choose(self):
        return random.choice(self.names)
    
    def alts_choose(self):
        return random.choice(self.alts)
        
    def bfactor_choose(self):
        return random.choice(self.bfactors)
        
    def occupancy_choose(self):
        return random.choice(self.occupancies)

    def generate_fragments(self):
        selection_fragments = []
        for seed in range(0,self.scale_n_tests_by):
            random.seed(seed)
            low,high = self.resseq_choose()
            fragments_comp = [
                f"resid {low} through {high}",
                f"resseq {low}:{high}",
                f"resname {self.resname_choose()}",
                f"chain {self.chain_choose()}",
                f"name {self.names_choose()}",
            ]
            fragments_other=[
                f"altloc '{self.alts_choose()}'",
                f"bfactor {self.operator_choose()} {self.bfactor_choose()}",
                f"occupancy {self.operator_choose()} {self.occupancy_choose()}",
            
            ]
            selection_fragments+=fragments_comp*3 + fragments_other # oversample compositional fragments
        self.selection_fragments = selection_fragments
        
    def filter_valid_fragments(self):
        # Test that the fragments are valid
        valid_fragments = []
        for fragment in self.selection_fragments:
            try:
                sel = self.model.selection(fragment).count(True)
                valid_fragments.append(fragment)
            except:
                self.log(f"Omitting invalid selection fragment: {fragment}")
        self.valid_fragments = valid_fragments

    def generate_random_expression(self,selection_fragments):
      
      # Randomly determine the number of fragments
      count = random.choice(self.expression_sizes)
      
      # Select fragments and operators
      choices = random.sample(selection_fragments, count)
      logic_choices = [random.choice(self.logic_tests) for _ in range(count - 1)] if count > 1 else []
      
      # Determine the number of parentheses groups
      n_parens = random.choice(range(count)) if count > 1 else 0
      
      # if count is 1, return the single fragment
      if count == 1:
        return choices[0]
    
      # if no parentheses are needed, just join fragments with logic operators
      if n_parens == 0:
        return " ".join(f"{frag} {logic}" for frag, logic in zip(choices[:-1], logic_choices)) + f" {choices[-1]}"
    
      # place parentheses and build the final expression
      expression = choices[0]
      open_parens = 0
    
      for i in range(1, count):
        if n_parens > 0 and random.choice([True, False]):
          expression = f"({expression}"
          open_parens += 1
          n_parens -= 1
        
        expression += f" {logic_choices[i-1]} {choices[i]}"
        
        if open_parens > 0 and random.choice([True, False]):
          expression += ")"
          open_parens -= 1
    
      # close parentheses
      expression += ")" * open_parens
    
      return expression
    def generate_expressions(self):
        expressions = []
        for seed in range(0,self.scale_n_tests_by*10):
            random.seed(seed)
            expression = self.generate_random_expression(self.valid_fragments)
            expressions.append(expression)
        self.selection_expressions = expressions

    def filter_expressions(self):
        # Filter expressions to be a specific number of atoms (not 0, not all)
        expressions_chosen = []
        for exp in self.selection_expressions:
            try:
                sel = self.model.selection(exp).count(True)
                if 5<sel<self.model.get_number_of_atoms():
                    expressions_chosen.append(exp)
            except:
                 self.log(f"Omitting invalid selection expression: {exp}")
        self.valid_expressions = expressions_chosen

    def test_selection_sequences(self):
        for seed in range(0,self.scale_n_tests_by):
            k = random.choice([2,3,4])
            selection = Selection(self.model)
            sample = random.sample(self.valid_expressions,k)
            for s in sample:
                selection = selection.select_from_string(s)
            selection.bool.count(True)
            try:
                selection.validate_debug()
            except:
                self.log("Failed selection sequence: ")
                for s in sample:
                    self.log(s)
                self.log()
                self.failures+=1
                
    def test_selection_reconstitution(self):
        for exp in self.valid_expressions:
            compatible, error, sel = PhenixParser.is_compatible_string(exp,self.model)
            if not compatible:
              self.log(f"Skipping: {exp} for reason: {error}")
            else:

              try:
                  phenix_string = self.reconstitute_string(sel)
                  sel_orig = self.model.selection(sel)
                  sel_reconstituted = self.model.selection(phenix_string)
                  assert (sel_orig==sel_reconstituted).count(True) == self.model.get_number_of_atoms(), (
                      "Input and output phenix strings do not select the same atoms")
              except:
                  
                  self.log("Failed, nonequivalent expressions detected:")
                  self.log(f"\tInitial expression: {exp}")
                  self.log(f"\tCompatible expression: {sel}")
                  self.log(f"\tReconstituted expression: {self.reconstitute_string(sel)}")
                  self.log()
                  #raise
                  self.failures+=1

    def reconstitute_string(self,string):
        parser = PhenixParser(string)
        parser.parse()
        phenix_string = parser.phenix_string
        return phenix_string
        
    def test_all(self):
        tester = self
        tester.generate_fragments()
        tester.filter_valid_fragments()
        tester.generate_expressions()
        tester.filter_expressions()
        #tester.test_selection_sequences()
        tester.test_selection_reconstitution()
        self.log(f"Tested {len(self.valid_expressions)} selection expressions.")
        self.log(f"Failures: {self.failures}")
        return self.failures==0



if __name__ == "__main__":
  app = QApplication(sys.argv)
  file = "/Users/user/Desktop/backup/3bet.cif"
  
  # Initialize your objects
  state = State.from_model_file(file)
  web_view = MolstarWebEngineView()
  web_page = WebEnginePage(web_view)
  web_view.setPage(web_page)
  controller = MolstarController(state, web_view=web_view)
  time.sleep(5)

  def run_tests_and_exit():
    tester = SelectionTester(model=state.first_model)  # Run tests
    tester.test_all()
    for exp in tester.valid_expressions:
      print(exp)
    # After tests are done, exit the application
    QApplication.quit()

  # Schedule the test to run once the event loop starts
  QTimer.singleShot(0, run_tests_and_exit)

  # Start the event loop
  sys.exit(app.exec_())

  # This will be printed after the application exits
  print("OK")