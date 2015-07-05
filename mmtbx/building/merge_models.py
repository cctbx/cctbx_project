from __future__ import division
import os,sys
import iotbx.pdb
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.maps.correlation
from scitbx.array_family import flex
from scitbx.matrix import col
from copy import deepcopy

# merge_models.py
# crossover models and pick best parts of each
#


master_phil = iotbx.phil.parse("""

  input_files {
    map_coeffs_file = None
      .type = path
      .help = File with map coefficients
      .short_caption = Map coefficients
      .style = bold file_type:hkl input_file process_hkl child:fobs:data_labels\
        child:space_group:space_group child:unit_cell:unit_cell

    map_coeffs_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of of map coefficients \
          to use
      .short_caption = Map coeffs label
      .style = bold renderer:draw_fobs_label_widget

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    pdb_in = None
      .type = path
      .help = Input PDB file to minimize
      .short_caption = Input PDB file

  }
  output_files {

    pdb_out = merged.pdb
      .type = path
      .help = Output PDB file (merged) 
      .short_caption = Output PDB file

  }
  crystal_info {
     resolution = None
       .type = float
       .help = High-resolution limit. Data will be truncated at this\
               resolution. If a map is supplied, it will be Fourier \
               filtered at this resolution. Required if input is a map and \
                only_original_map is not set.
       .short_caption = High-resolution limit
       .style = resolution
     space_group = None
       .type = space_group
       .short_caption = Space Group
       .help = Space group (normally read from the data file)
     unit_cell = None
       .type = unit_cell
       .short_caption = Unit Cell
       .help = Unit Cell (normally read from the data file)

     scattering_table = wk1995  it1992  *n_gaussian  neutron electron
       .type = choice
       .help = Scattering table to use
       .short_caption = Scattering table
  }
  crossover {

     minimum_length = 2
       .type = int
       .short_caption = Minimum length of a crossover
       .help = Minimum length of a crossover

     maximum_fraction = 0.5
       .type = float
       .short_caption = Maximum replacement fraction
       .help = Maximum replacement fraction

     dist_max = 1.0
       .type = float
       .short_caption = Crossover distance
       .help = Maximum distance between atoms where crossover is to occur

     minimum_matching_atoms = 3
       .type = int
       .short_caption = Minimum number of matching atoms to cross over 
       .help =  Minimum number of matching atoms to cross over

     crossover_atom =  CA
       .type = str 
       .short_caption = Crossover atom
       .help = Atom where crossovers will occur

     smoothing_window = 5
       .type = int 
       .short_caption = Smoothing window 
       .help = Smoothing window. The residue CC values will be smoothed with \
              this window in calculating the optimal crossover 

  }
  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output
  }
""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):
  command_line = iotbx.phil.process_command_line_with_files(
    reflection_file_def="input_files.map_coeffs_file",
    map_file_def="input_files.map_file",
    pdb_file_def="input_files.pdb_in",
    args=args,
    master_phil=master_phil)
  params = command_line.work.extract()
  print >>out,"\nMerge_models: Take parts of multiple models to construct one model\n"
  master_phil.format(python_object=params).show(out=out)
  return params

class model_object:
  def __init__(self,source_id=None,source_list=None,cc_dict=None,
      crossover_dict=None,
      minimum_length=1,
      maximum_fraction=0.5):
    if source_id is not None:
       self.source_list=cc_dict[source_id].size()*[source_id]
    else:
      self.source_list=source_list # list of models as source for each residue
    assert self.source_list is not None

    self.size=len(self.source_list)
    self.cc_dict=cc_dict
    self.minimum_length=minimum_length
    self.maximum_fraction=maximum_fraction
    self.crossover_dict=crossover_dict
    self.reset_score()

  def show_summary(self,out=sys.stdout):
    print >>out,"\nModel with %d sites and score of %7.2f" %(
     len(self.source_list),self.score)
    print >>out," ".join(self.source_list).replace("  "," ")

  def is_allowed_crossover(self,i=None,other=None):
    # return True if a crossover at position i to model other is allowed
    if other.source_list[i] in \
        self.crossover_dict.get(i,{}).get(self.source_list[i],[]):
      return True
    return False

  def customized_copy(self):
    new_model=model_object(
     source_list=deepcopy(self.source_list),
     minimum_length=self.minimum_length,
     maximum_fraction=self.maximum_fraction,
     cc_dict=self.cc_dict,
     crossover_dict=self.crossover_dict)
    return new_model

  def optimize_with_others(self,others=None):
    found=True
    best_model=self
    cycle=0
    while found:
      cycle+=1
      found=False
      for other_model in others:
          new_model=best_model.select_best_from_other(other_model)
          if not new_model: continue
          if best_model is None or new_model.get_score()>best_model.get_score():
            best_model=new_model
            found=True
    return best_model

  def select_best_from_other(self,other=None):
    # find the part of other that can best replace part of this one using
    # allowed crossovers
    best_score=None
    best_model=None
    for i1 in xrange(self.size):
      if not self.is_allowed_crossover(i1,other): continue
      for i2 in xrange(i1+self.minimum_length+1,self.size):
        if not self.is_allowed_crossover(i2,other): continue
        if float(i2-i1+1)/self.size > self.maximum_fraction: break
          # test replacing self with i1 to i2 of other
        test_model=self.customized_copy()
        for i in xrange(i1,i2+1):
            test_model.source_list[i]=other.source_list[i]
        test_model.reset_score()
        if best_score is None or test_model.get_score()>best_score:
          best_score=test_model.get_score()
          best_model=test_model
    return best_model
         
            
  def reset_score(self):
    self.score=None
    self.get_score()

  def get_score(self):
    if self.score is not None:
       return self.score

    # return None if any stretch of residues from 1 model is shorter than
    #  self.minimum_length
    last_id=None
    n=0
    for i in xrange(self.size):
      if last_id is None: last_id=self.source_list[i]
      if self.source_list[i]==last_id:
        n+=1
      elif n < self.minimum_length:
        return None
      else:
        n=1
        last_id=self.source_list[i]
    if n>0 and n<self.minimum_length:
      return None


    
    # sum up CC values at each residue
    score=0.
    for i in xrange(self.size):
      score+=self.cc_dict[self.source_list[i]][i]
    self.score=score
    return score

def get_crossover_dict(
      n_residues=None,
      hierarchy=None,chain_id=None,
      crossover_atom=None,
      minimum_matching_atoms=None,
      dist_max=None,
      verbose=None,
      out=sys.stdout):
  crossover_dict={}  # Allowed crossover for [position][id1][id2]
  print >>out, "\nMaking a table of allowed crossovers"

  # select out just the crossover atoms...
  atom_selection="name %s" %(crossover_atom)
  asc=hierarchy.atom_selection_cache()
  sel=asc.selection(string = atom_selection)
  sel_hierarchy=hierarchy.select(sel)

  dist_max_sq=dist_max**2
  used_model_ids=[]
  for model1 in sel_hierarchy.models():
    used_model_ids.append(model1.id)
    for chain1 in model1.chains():
      if chain1.id != chain_id: continue
      xyz1=chain1.atoms().extract_xyz()
      for model2 in sel_hierarchy.models():
        if model2.id in used_model_ids: continue # already did it
        for chain2 in model2.chains():
          if chain2.id != chain_id: continue
          xyz2=chain2.atoms().extract_xyz()
          for i in xrange(xyz1.size()):
            x1=col(xyz1[i])
            x2=col(xyz2[i])
            dd=(x1-x2).norm_sq()
            if dd<= dist_max_sq:  # can crossover here
              if not i in crossover_dict.keys(): crossover_dict[i]={}
              if not model1.id in crossover_dict[i].keys():
                crossover_dict[i][model1.id]=[]
              if not model2.id in crossover_dict[i].keys():
                crossover_dict[i][model2.id]=[]
              if not model2.id in crossover_dict[i][model1.id]:
                  crossover_dict[i][model1.id].append(model2.id)
              if not model1.id in crossover_dict[i][model2.id]:
                  crossover_dict[i][model2.id].append(model1.id)
            
  # Now remove where the number
  #  of crossover atoms in a row that match is less than minimum_matching_atoms

  if minimum_matching_atoms > 2:
    offset_n=minimum_matching_atoms//2
    offset_range=[]
    for n in xrange(-offset_n,offset_n+1):
      if n != 0: offset_range.append(n)
    delete_dict={}
    for i in xrange(n_residues):
      if not i in crossover_dict.keys(): continue
      for id1 in crossover_dict[i].keys():
        for id2 in crossover_dict[i][id1]:
          # check to see if i-1 and i+1 are both ok (if not off the ends)
          for offset in offset_range:
            i1=min(n_residues-1,max(0,i+offset))
            if not id2 in crossover_dict.get(i1,{}).get(id1,[]):
              if not i in delete_dict.keys(): delete_dict[i]={}
              if not id1 in delete_dict[i].keys(): delete_dict[i][id1]=[]
              if not id2 in delete_dict[i][id1]:delete_dict[i][id1].append(id2)

    for i in xrange(n_residues):
      if not i in crossover_dict.keys(): continue
      for id1 in crossover_dict[i].keys():
        new_list=[]
        for id2 in crossover_dict[i][id1]:
          if not id2 in delete_dict.get(i,{}).get(id1,[]):
            new_list.append(id2)
        crossover_dict[i][id1]=new_list

  # Now add all ends to crossover (always ok) 

  for pos in [0,n_residues-1]:
    if not pos in crossover_dict.keys():
      crossover_dict[pos]={}
    for id1 in used_model_ids:
      for id2 in used_model_ids:
        if id1==id2: continue
        if not id1 in crossover_dict[pos].keys(): 
          crossover_dict[pos][id1]=[]
        if not id2 in crossover_dict[pos][id1]: 
          crossover_dict[pos][id1].append(id2)
 
   
  if verbose: 
    i_list=crossover_dict.keys()
    i_list.sort()
    for i in i_list:
      print >>out,"\nAllowed crossovers at position %d" %(i)
      id_list=crossover_dict[i].keys()
      id_list.sort()
      print "Crossover pairs:",
      for id in id_list:
        second_id_list=crossover_dict[i][id]
        second_id_list.sort()
        for second_id in second_id_list:
          print "%s-%s" %(id,second_id),
      print

  return crossover_dict


def get_cc_dict(hierarchy=None,chain_id=None,cc_calculator=None,out=sys.stdout):
  cc_dict={}
  print >>out, "\nMaking a table of residue CC values for chain %s" %(chain_id)
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id != chain_id: continue
      #print "\nMODEL %s CHAIN %s:\n" %(model.id,chain.id)
      cc_list=flex.double()
      cc_dict[model.id]=cc_list
      for rg in chain.residue_groups():
        cc = cc_calculator.cc(selection=rg.atoms().extract_i_seq())
        #print >> out, "  chain id: %s resid %s: %6.4f"%(
            #rg.parent().id, rg.resid(), cc)
        cc_list.append(cc)

  # check to make sure all are same
  std_size=None
  for model_id in cc_dict.keys():
    if std_size is None: std_size=cc_dict[model_id].size()
    if cc_dict[model_id].size()!=std_size:
      raise Sorry(
       "All copies of each chain must have the same size (%d != %d)" %(
         std_size,cc_dict[model_id].size()))
  return cc_dict

def smooth_cc_values(cc_dict=None,
       smoothing_window=None,verbose=None,out=sys.stdout):
  smoothed_cc_dict={}
  delta=(smoothing_window+1)//2
  for id in cc_dict.keys():
    cc_list=cc_dict[id]
    smoothed_cc_list=flex.double()
    for i in xrange(cc_list.size()):
      r=cc_list[max(0,i-delta):min(cc_list.size(),i+delta+1)]
      smoothed_cc_list.append(r.min_max_mean().mean)
    smoothed_cc_dict[id]=smoothed_cc_list

  keys=smoothed_cc_dict.keys()
  keys.sort()
  if verbose:
      for key in keys:
        print >>out,"ID:  %s " %(key)
        print >>out,"Position   Unsmoothed  Smoothed"
        for i in xrange(cc_dict[key].size()):
         print >>out,"  %d     %7.2f     %7.2f" %(
           i,cc_dict[key][i],smoothed_cc_dict[key][i])

  return smoothed_cc_dict

def remove_ter(text): # remove blank lines and TER records
  new_lines=[]
  for line in flex.split_lines(text):
    if not line.replace(" ",""): continue
    if line.startswith("TER"): continue
    new_lines.append(line)
  return "\n".join(new_lines)

def run(args,
    map_data=None,
    map_coeffs=None,
    pdb_inp=None,
    states=None,
    pdb_string=None,
    crystal_symmetry=None,
    out=sys.stdout):

  # Get the parameters
  params=get_params(args=args,out=out)

  # Get map_data if not present
  if not map_data:
    if not map_coeffs:
      from mmtbx.building.minimize_chain import get_map_coeffs
      map_coeffs=get_map_coeffs(
        map_coeffs_file=params.input_files.map_coeffs_file,
        map_coeffs_labels=params.input_files.map_coeffs_labels)
    if not map_coeffs:
      raise Sorry("Need map_coeffs_file")

    fft_map = map_coeffs.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()

  map_data=map_data.as_double()

  if map_coeffs and not crystal_symmetry:
    crystal_symmetry=map_coeffs.crystal_symmetry()

  if map_coeffs and not params.crystal_info.resolution:
    params.crystal_info.resolution=map_coeffs.d_min()

  assert crystal_symmetry is not None

  # Get the starting model
  if pdb_inp is None:
    if not pdb_string:
      if params.input_files.pdb_in:
        print >>out,"Taking models from %s" %(params.input_files.pdb_in)
        pdb_string=open(params.input_files.pdb_in).read()
      elif states:
        print >>out,"Taking models from input states"
        pdb_string=states.root.as_pdb_string()
      else:
        raise Sorry("Need an input PDB file")
    pdb_inp=iotbx.pdb.input(source_info=None, lines = pdb_string)
    if not pdb_inp.crystal_symmetry(): # get it
      cryst1_line=iotbx.pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry)
      from cStringIO import StringIO
      f=StringIO()
      print >>f, cryst1_line
      print >>f,pdb_string
      pdb_string=f.getvalue()
      pdb_inp=iotbx.pdb.input(source_info=None, lines = pdb_string)

  if not pdb_inp:
    raise Sorry("Need a model or models")

  hierarchy = pdb_inp.construct_hierarchy()
  n_models=0
  for model in hierarchy.models():
    n_models+=1

  if n_models==1:  # nothing to do
    return hierarchy

  xrs = pdb_inp.xray_structure_simple(crystal_symmetry=crystal_symmetry)
  xrs.scattering_type_registry(table = params.crystal_info.scattering_table)

  if not params.crystal_info.resolution:
    from cctbx import maptbx
    params.crystal_info.resolution=maptbx.resolution_from_map_and_model(
      map_data=map_data, xray_structure=xrs)

  print >>out,"\nResolution limit: %7.2f" %(params.crystal_info.resolution)
  print >>out,"\nSummary of input models"
  xrs.show_summary(f=out, prefix="  ")

  print >>out, "\nReady with %d models and map" %(n_models)
  # Get CC by residue for each model and map

  cc_calculator = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
    xray_structure = xrs,
    map_data       = map_data,
    d_min          = params.crystal_info.resolution)



  # Run through chains separately
  chain_id_list=[]
  for model in hierarchy.models():
    for chain in model.chains():
      if not chain.id in chain_id_list: chain_id_list.append(chain.id)

  for chain_id in chain_id_list:



    # get CC values for all residues
    cc_dict=get_cc_dict(hierarchy=hierarchy,chain_id=chain_id,
     cc_calculator=cc_calculator,out=out)


    # smooth CC values with window of smoothing_window
    smoothed_cc_dict=smooth_cc_values(cc_dict=cc_dict,
       smoothing_window=params.crossover.smoothing_window,
       verbose=params.control.verbose,out=out)

    # figure out all the places where crossover can occur.
    # this version: only crossover where corresponding residues are within
    # dist_max of each other.
    n_residues=cc_dict[cc_dict.keys()[0]].size()

    crossover_dict=get_crossover_dict(
      n_residues=n_residues,
      hierarchy=hierarchy,chain_id=chain_id,
      crossover_atom=params.crossover.crossover_atom, 
      dist_max=params.crossover.dist_max,
      minimum_matching_atoms=params.crossover.minimum_matching_atoms,
      verbose=params.control.verbose,out=out)


    # Now we are ready to identify the best composite model...
    # A composite has reside 0 from model x, residue 1 from model y etc.
    # Each change from model a to model b between residues i and i+1 must have
    #  a crossover between a and b at either residue i or i+1 

    keys=cc_dict.keys()
    keys.sort()

    working_model_list=[]
    for key in keys:
      working_model=model_object(source_id=key,cc_dict=cc_dict,
         crossover_dict=crossover_dict,
         minimum_length=params.crossover.minimum_length,
         maximum_fraction=params.crossover.maximum_fraction)
      if params.control.verbose:
        working_model.show_summary(out=out)
      working_model_list.append(working_model)


    # Go through all the working models and cross them with other models to
    #  optimize...Then take all the best and cross...

    best_model=working_model_list[0]
    found=True
    while found:
      found=False
      new_working_model_list=[]
      new_best=best_model
      for working_model in working_model_list:
        others=[]
        for m in working_model_list:
          if not working_model==m:  others.append(m)
        new_working_model=working_model.optimize_with_others(others=others)
        new_working_model_list.append(new_working_model)
        if new_working_model.get_score()>best_model.get_score():
          new_best=new_working_model
      if new_best.get_score()>best_model.get_score():
        if params.control.verbose:
          print "NEW BEST SCORE: %7.2f" %(new_best.get_score())
          new_best.show_summary(out=out)
        best_model=new_best
        found=True
        


    # Write out best_model

    # note residue values
    for model in hierarchy.models():
      for chain in model.chains():
        if chain.id != chain_id: continue
        residue_list=[]
        for rg in chain.residue_groups():
          residue_list.append(rg.resseq)
    residue_list.sort()
    assert len(best_model.source_list)==len(residue_list)

    from cStringIO import StringIO
    f=StringIO()
    for i in xrange(len(residue_list)): #
      atom_selection="model %s and resseq %s" %(best_model.source_list[i],
          residue_list[i])
      asc=hierarchy.atom_selection_cache()
      sel=asc.selection(string = atom_selection)
      sel_hierarchy=hierarchy.select(sel)
      print >>f,remove_ter(sel_hierarchy.as_pdb_string())
    pdb_string=f.getvalue()
    pdb_inp=iotbx.pdb.input(source_info=None, lines = pdb_string)
    xray_structure=pdb_inp.xray_structure_simple()
    pdb_hierarchy=pdb_inp.construct_hierarchy()

    if params.output_files.pdb_out:
      f=open(params.output_files.pdb_out,'w')
      print >>f,pdb_string
      print "Final model is in: \n%s" %(f.name)
      f.close()

    return pdb_hierarchy, xray_structure

if   (__name__ == "__main__"):
  args=sys.argv[1:]
  run(args=args)
