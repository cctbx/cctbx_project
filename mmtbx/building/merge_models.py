from __future__ import absolute_import, division, print_function
import sys,os
import iotbx.pdb
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.maps.correlation
from scitbx.array_family import flex
from scitbx.matrix import col
from copy import deepcopy
from libtbx import adopt_init_args
from cctbx.maptbx import resolution_from_map_and_model
from six.moves import range
from six.moves import cStringIO as StringIO

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

    pdb_in_file = None
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

     max_keep = 10
       .type = int
       .short_caption = Max keep
       .help = Max keep. Number of solutions to carry along during optimization

     max_regions_to_test = 10
       .type = int
       .short_caption = Max regions to test
       .help = Maximum number of regions within a chain to test for crossover \
               with another chain.  Sorted on smoothed differences in local CC

     max_ends_per_region = 5
       .type = int
       .short_caption = Max ends per region
       .help = Maximum number of ends (left and right) to test for each \
               potential crossover region.

     minimum_improvement = 0.01
       .type = float
       .short_caption = Minimum improvement
       .help = Minimum improvement to keep a crossover

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
    pdb_file_def="input_files.pdb_in_file",
    args=args,
    master_phil=master_phil)
  params = command_line.work.extract()
  print("\nMerge_models: Take parts of multiple models to construct one model\n", file=out)
  master_phil.format(python_object=params).show(out=out)
  return params

class model_object:
  def __init__(self,
      source_id=None,
      source_list=None,
      cc_dict=None,
      smoothed_cc_dict=None,
      crossover_dict=None,
      minimum_length=0,
      minimum_improvement=0.01,
      max_regions_to_test=15,
      max_ends_per_region=5,
      maximum_fraction=0.5):

    adopt_init_args(self, locals())

    if self.source_id is not None: # list of models as source for each residue
       self.source_list=cc_dict[self.source_id].size()*[self.source_id]
    assert self.source_list is not None

    self.size=len(self.source_list)
    self.reset_score()

  def show_summary(self,out=sys.stdout):
    print("\nModel with %d sites and score of %7.2f" %(
     len(self.source_list),self.score), file=out)
    print(" ".join(self.source_list).replace("  "," "), file=out)

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
     minimum_improvement=self.minimum_improvement,
     maximum_fraction=self.maximum_fraction,
     max_regions_to_test=self.max_regions_to_test,
     max_ends_per_region=self.max_ends_per_region,
     cc_dict=self.cc_dict,
     smoothed_cc_dict=self.smoothed_cc_dict,
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
          new_model=best_model.select_best_from_other_smooth(other_model)
          if not new_model: continue
          if best_model is None or new_model.get_score()>best_model.get_score():
            best_model=new_model
            found=True
    return best_model

  def select_best_from_other_smooth(self,
    other=None):

    # Use smoothed cc_dict to identify location where
    # biggest difference can be found and focus on that

    difference_list=[]
    for i1 in range(self.size):
      source_self=self.source_list[i1]
      source_other=other.source_list[i1]
      difference_list.append( [
       self.cc_dict[source_other][i1]-self.cc_dict[source_self][i1],
       i1])

    difference_list.sort()
    difference_list.reverse()

    best_score=self.get_score()
    original_score=best_score
    best_model=self

    # Now work down difference list until we get something useful.
    # For each one, find how far in either direction we can to go
    #  maximize the score after crossover

    for dd,i in difference_list[:self.max_regions_to_test]:
      if dd<self.minimum_improvement: continue

      #Now figure out where we can cross over on either side of i
      allowed_left_crossovers=[]
      allowed_right_crossovers=[]
      for ib in range(self.size):
        if not self.is_allowed_crossover(i,other): continue
        if ib <i: allowed_left_crossovers.append(ib)
        if ib >i: allowed_right_crossovers.append(ib)

      if not allowed_left_crossovers or not allowed_right_crossovers: continue

      # find best to left and to right, up to max_ends_per_region

      # if possible, find a right-crossover that is minimum_length from center
      i2=allowed_right_crossovers[0]
      for test_i2 in allowed_right_crossovers:
        if test_i2-i>self.minimum_length: # take it; it is long enough
          i2=test_i2
          break

      # find best start (i1, on left)
      best_i=None
      best_i_score=None
      for i1 in allowed_left_crossovers[-self.max_ends_per_region:]:
        if float(i2-i1+1)/self.size > self.maximum_fraction: continue
        # test replacing self with i1 to i2 of other
        test_model=self.customized_copy()
        for i in range(i1,i2+1):
          test_model.source_list[i]=other.source_list[i]
        test_model.reset_score()
        if best_i_score is None or \
           (test_model.get_score() and test_model.get_score()>best_i_score):
          best_i_score=test_model.get_score()
          best_i=i1
      i1=best_i
      if i1 is None:
        continue

      # Now find best end (right side ;i2)
      best_i=None
      best_i_score=None

      for i2 in allowed_right_crossovers[:self.max_ends_per_region]:
        if float(i2-i1+1)/self.size > self.maximum_fraction: continue
        # test replacing self with i1 to i2 of other
        test_model=self.customized_copy()
        for i in range(i1,i2+1):
          test_model.source_list[i]=other.source_list[i]
        test_model.reset_score()
        if best_i_score is None or \
           (test_model.get_score() and test_model.get_score()>best_i_score):
          best_i_score=test_model.get_score()
          best_i=i2

          #  save if best overall
          if test_model.get_score() and \
             test_model.get_score()>best_score+self.minimum_improvement:
            best_score=test_model.get_score()
            best_model=test_model
      if best_model.get_score()>original_score+self.minimum_improvement:
        return best_model  # take it (and skip other possibilities)


  def reset_score(self):
    self.score=None
    self.get_score()

  def get_score(self):
    if self.score is not None:
       return self.score

    # return -999 if any stretch of residues from 1 model is shorter than
    #  self.minimum_length
    last_id=None
    n=0
    for i in range(self.size):
      if last_id is None: last_id=self.source_list[i]
      if self.source_list[i]==last_id:
        n+=1
      elif n < self.minimum_length:
        return -999.
      else:
        n=1
        last_id=self.source_list[i]
    if n>0 and n<self.minimum_length:
      return -999.



    # sum up CC values at each residue
    score=0.
    for i in range(self.size):
      score+=self.cc_dict[self.source_list[i]][i]
    self.score=score
    return score

def get_atom_selection(chain_id=None,model_id=None,resseq_sel=None,
   start_resno=None,end_resno=None):

  if chain_id and chain_id.replace(" ",""):
    chain_sel=" chain %s " %(chain_id)
  else:
    chain_sel=""
  if model_id and model_id.replace(" ",""):
    model_sel="model %s " %(model_id)
  else:
    model_sel=""
  if chain_sel and model_sel:
    and_sel=" and "
  else:
    and_sel=""
  atom_selection="%s %s %s" %(model_sel,and_sel,chain_sel)

  if resseq_sel and resseq_sel.replace(" ",""):
    if atom_selection.replace(" ",""):
      atom_selection="%s and resseq %s" %(atom_selection,resseq_sel)
    else:
      atom_selection="resseq %s" %(resseq_sel)

  elif start_resno is not None and end_resno is not None:
    if atom_selection.replace(" ",""):
      atom_selection="%s and resseq %d:%d" %(atom_selection,start_resno,end_resno)
    else:
      atom_selection="resseq %d:%d" %(start_resno,end_resno)

  return atom_selection

def get_crossover_dict(
      n_residues=None,
      hierarchy=None,
      crossover_atom=None,
      minimum_matching_atoms=None,
      dist_max=None,
      verbose=None,
      out=sys.stdout):
  crossover_dict={}  # Allowed crossover for [position][id1][id2]
  print("\nMaking a table of allowed crossovers", file=out)

  # select out just the crossover atoms...

  atom_selection="name %s " %(crossover_atom)

  asc=hierarchy.atom_selection_cache()
  sel=asc.selection(string = atom_selection)
  sel_hierarchy=hierarchy.select(sel)

  dist_max_sq=dist_max**2
  used_model_ids=[]
  for model1 in sel_hierarchy.models():
    used_model_ids.append(model1.id)
    for chain1 in model1.chains():
      xyz1=chain1.atoms().extract_xyz()
      for model2 in sel_hierarchy.models():
        if model2.id in used_model_ids: continue # already did it
        for chain2 in model2.chains():
          xyz2=chain2.atoms().extract_xyz()
          if xyz1.size()!=xyz2.size():
            print("\nSize of chain " +\
            "'%s' model '%s' (%d) is different than chain '%s' model '%s' (%d) " %(
              chain1.id,model1.id,xyz1.size(),chain2.id,model2.id,xyz2.size()), file=out)
            assert xyz1.size()==xyz2.size()

          for i in range(xyz1.size()):
            x1=col(xyz1[i])
            x2=col(xyz2[i])
            dd=(x1-x2).norm_sq()
            if dd<= dist_max_sq:  # can crossover here
              if not i in crossover_dict: crossover_dict[i]={}
              if not model1.id in crossover_dict[i]:
                crossover_dict[i][model1.id]=[]
              if not model2.id in crossover_dict[i]:
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
    for n in range(-offset_n,offset_n+1):
      if n != 0: offset_range.append(n)
    delete_dict={}
    for i in range(n_residues):
      if not i in crossover_dict: continue
      for id1 in crossover_dict[i]:
        for id2 in crossover_dict[i][id1]:
          # check to see if i-1 and i+1 are both ok (if not off the ends)
          for offset in offset_range:
            i1=min(n_residues-1,max(0,i+offset))
            if not id2 in crossover_dict.get(i1,{}).get(id1,[]):
              if not i in delete_dict: delete_dict[i]={}
              if not id1 in delete_dict[i]: delete_dict[i][id1]=[]
              if not id2 in delete_dict[i][id1]:delete_dict[i][id1].append(id2)

    for i in range(n_residues):
      if not i in crossover_dict: continue
      for id1 in crossover_dict[i]:
        new_list=[]
        for id2 in crossover_dict[i][id1]:
          if not id2 in delete_dict.get(i,{}).get(id1,[]):
            new_list.append(id2)
        crossover_dict[i][id1]=new_list

  # Now add all ends to crossover (always ok)

  for pos in [0,n_residues-1]:
    if not pos in crossover_dict:
      crossover_dict[pos]={}
    for id1 in used_model_ids:
      for id2 in used_model_ids:
        if id1==id2: continue
        if not id1 in crossover_dict[pos]:
          crossover_dict[pos][id1]=[]
        if not id2 in crossover_dict[pos][id1]:
          crossover_dict[pos][id1].append(id2)


  if verbose:
    i_list=list(crossover_dict.keys())
    i_list.sort()
    for i in i_list:
      print("\nAllowed crossovers at position %d" %(i), file=out)
      id_list=list(crossover_dict[i].keys())
      id_list.sort()
      print("Crossover pairs:", end=' ')
      for id in id_list:
        second_id_list=crossover_dict[i][id]
        second_id_list.sort()
        for second_id in second_id_list:
          print("%s-%s" %(id,second_id), end=' ')
      print()

  return crossover_dict


def get_cc_dict(hierarchy=None,crystal_symmetry=None,
  map_data=None,d_min=None,
  table=None,out=sys.stdout):

  cc_dict={}
  print("\nMaking a table of residue CC values", file=out)
  cryst1_line=iotbx.pdb.format_cryst1_record(crystal_symmetry=crystal_symmetry)

  # select the model and chain we are interested in
  for model in hierarchy.models():

    f=StringIO()
    atom_selection=get_atom_selection(model_id=model.id)
    asc=hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection)
    sel_hierarchy=hierarchy.select(sel)
    sel_hierarchy.atoms().reset_i_seq()
    xrs = sel_hierarchy.extract_xray_structure(
      crystal_symmetry=crystal_symmetry)
    xrs.scattering_type_registry(table = table)

    cc_calculator=mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      xray_structure = xrs,
      map_data       = map_data,
      d_min          = d_min)

    for m in sel_hierarchy.models():
      for chain in m.chains():
        cc_list=flex.double()
        cc_dict[model.id]=cc_list
        for rg in chain.residue_groups():
          cc = cc_calculator.cc(selection=rg.atoms().extract_i_seq())
          #print >> out, "  chain id: %s resid %s: %6.4f"%( rg.parent().id, rg.resid(), cc)
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
    for i in range(cc_list.size()):
      r=cc_list[max(0,i-delta):min(cc_list.size(),i+delta+1)]
      smoothed_cc_list.append(r.min_max_mean().mean)
    smoothed_cc_dict[id]=smoothed_cc_list

  keys=list(smoothed_cc_dict.keys())
  keys.sort()
  if verbose:
      for key in keys:
        print("ID:  %s " %(key), file=out)
        print("Position   Unsmoothed  Smoothed", file=out)
        for i in range(cc_dict[key].size()):
         print("  %d     %7.2f     %7.2f" %(
           i,cc_dict[key][i],smoothed_cc_dict[key][i]), file=out)

  return smoothed_cc_dict

# NOTE: Match defaults here and in params at top of file
#     : copy from defaults if params is not None below
#     : See explanations of parameters in params at top of file
def run(
    params=None, # params for running from command line
    map_data=None,  # map_data, as_double()
    pdb_inp=None,
    pdb_hierarchy=None,
    crystal_symmetry=None,
    resolution=None,
    scattering_table='n_gaussian',
    smoothing_window=5,
    crossover_atom='CA',
    minimum_matching_atoms=3,
    minimum_length=2,
    dist_max=1.0,
    minimum_improvement=0.01,
    max_regions_to_test=10,
    max_ends_per_region=5,
    maximum_fraction=0.5,
    max_keep=10,
    map_coeffs_file=None,map_coeffs_labels=None,
    pdb_in_file=None,
    pdb_out=None,
    verbose=None,
    out=sys.stdout):

  if out is None: out=sys.stdout # explode and refine calls it this way

  # get info from params if present
  if params:
     verbose=params.control.verbose
     map_coeffs_file=params.input_files.map_coeffs_file
     map_coeffs_labels=params.input_files.map_coeffs_labels
     pdb_in_file=params.input_files.pdb_in_file
     resolution=params.crystal_info.resolution
     scattering_table=params.crystal_info.scattering_table
     smoothing_window=params.crossover.smoothing_window
     crossover_atom=params.crossover.crossover_atom
     minimum_matching_atoms=params.crossover.minimum_matching_atoms
     minimum_length=params.crossover.minimum_length
     dist_max=params.crossover.dist_max
     minimum_improvement=params.crossover.minimum_improvement
     max_regions_to_test=params.crossover.max_regions_to_test
     max_ends_per_region=params.crossover.max_ends_per_region
     maximum_fraction=params.crossover.maximum_fraction
     max_keep=params.crossover.max_keep
     pdb_out=params.output_files.pdb_out

  # Consistency checks
  if(pdb_hierarchy is not None):
    assert pdb_in_file is None
    assert pdb_inp is None
    assert crystal_symmetry is not None
    # XXX more checks here!

  # Get map_data if not present
  if not map_data:
    if not map_coeffs_file or not os.path.isfile(map_coeffs_file):
      raise Sorry("Cannot find the map_coeffs_file '%s'" %(
        str(map_coeffs_file)))
    from mmtbx.building.minimize_chain import get_map_coeffs
    map_coeffs=get_map_coeffs(map_coeffs_file,
        map_coeffs_labels=map_coeffs_labels)

    fft_map = map_coeffs.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    map_data=map_data.as_double()
    if map_coeffs and not crystal_symmetry:
      crystal_symmetry=map_coeffs.crystal_symmetry()
    if map_coeffs and not resolution:
      resolution=map_coeffs.d_min()

  # Get the starting model
  if(pdb_hierarchy is None):
    if pdb_inp is None:
      if not pdb_in_file or not os.path.isfile(pdb_in_file):
        raise Sorry("Cannot read input PDB file '%s'" %(
          str(pdb_in_file)))
      else:
        print("Taking models from %s" %(pdb_in_file), file=out)
        pdb_string=open(pdb_in_file).read()
      pdb_inp=iotbx.pdb.input(source_info=None, lines = pdb_string)
      if pdb_inp is None:
        raise Sorry("Need a model or models")
    if not crystal_symmetry:
      crystal_symmetry=pdb_inp.crystal_symmetry()
    assert crystal_symmetry is not None
    hierarchy = pdb_inp.construct_hierarchy()
  else:
    hierarchy = pdb_hierarchy # XXX FIXME
  n_models=0
  for model in hierarchy.models():
    n_models+=1

  if n_models==1:  # nothing to do
    return hierarchy

  xrs = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  xrs.scattering_type_registry(table=scattering_table)
  if not resolution:
    from cctbx import maptbx
    resolution=maptbx.resolution_from_map_and_model.run(
      map_data=map_data, xray_structure=xrs).d_min
  if(resolution is None):
    raise Sorry("Resolution is required")
  print("\nResolution limit: %7.2f" %(resolution), file=out)
  print("\nSummary of input models", file=out)
  xrs.show_summary(f=out, prefix="  ")

  print("\nReady with %d models and map" %(n_models), file=out)
  # Get CC by residue for each model and map

  chain_id_and_resseq_list=[] # Instead set up chain_id and resseq (range)
  from mmtbx.secondary_structure.find_ss_from_ca import \
      split_model
  model_list=split_model(hierarchy=hierarchy,only_first_model=True)
  for m in model_list:
    h=m.hierarchy
    first_resno=h.first_resseq_as_int()
    last_resno=h.last_resseq_as_int()
    chain_id=h.first_chain_id()
    residue_range=[first_resno,last_resno]
    chain_id_and_resseq=[chain_id,residue_range]
    if not chain_id_and_resseq in chain_id_and_resseq_list:
       chain_id_and_resseq_list.append(chain_id_and_resseq)

  # Run through chains separately
  # NOTE: All models of each chain must match exactly

  # Save composite model, chain by chain
  composite_model_stream=StringIO()
  sel_ph_list = []

  for chain_id_and_resseq in chain_id_and_resseq_list:
    f=StringIO()
    chain_id,[start_resno,end_resno]=chain_id_and_resseq
    atom_selection=get_atom_selection(chain_id=chain_id,
      start_resno=start_resno,end_resno=end_resno)
    asc=hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection)
    sel_hierarchy=hierarchy.select(sel)
    sel_hierarchy.atoms().reset_i_seq()
    ph=sel_hierarchy

    print("\nWorking on chain_id='%s' resseq %d:%d\n" %(
       chain_id_and_resseq[0],chain_id_and_resseq[1][0],chain_id_and_resseq[1][1]), file=out)

    # get CC values for all residues
    cc_dict=get_cc_dict(hierarchy=ph,map_data=map_data,d_min=resolution,
     crystal_symmetry=crystal_symmetry,
     table=scattering_table,out=out)

    # smooth CC values with window of smoothing_window
    smoothed_cc_dict=smooth_cc_values(cc_dict=cc_dict,
       smoothing_window=smoothing_window,
       verbose=verbose,out=out)

    # figure out all the places where crossover can occur.
    # FIXME: order of keys changes in py2/3 vthis could be bad. No all are same.
    n_residues=cc_dict[list(cc_dict.keys())[0]].size()

    crossover_dict=get_crossover_dict(
      n_residues=n_residues,
      hierarchy=ph,
      crossover_atom=crossover_atom,
      dist_max=dist_max,
      minimum_matching_atoms=minimum_matching_atoms,
      verbose=verbose,out=out)

    # Now we are ready to identify the best composite model...
    # A composite has reside 0 from model x, residue 1 from model y etc.
    # Each change from model a to model b between residues i and i+1 must have
    #  a crossover between a and b at either residue i or i+1

    keys=list(cc_dict.keys())
    keys.sort()

    sorted_working_model_list=[]
    for key in keys:
      working_model=model_object(source_id=key,
         cc_dict=cc_dict,
         smoothed_cc_dict=smoothed_cc_dict,
         crossover_dict=crossover_dict,
         minimum_length=minimum_length,
         minimum_improvement=minimum_improvement,
         max_regions_to_test=max_regions_to_test,
         max_ends_per_region=max_ends_per_region,
         maximum_fraction=maximum_fraction)
      if verbose:
        working_model.show_summary(out=out)
      sorted_working_model_list.append(
        [working_model.get_score(),working_model])
    sorted_working_model_list.sort()
    sorted_working_model_list.reverse()
    sorted_working_model_list=\
       sorted_working_model_list[:max_keep]
    working_model_list=[]
    for s,m in sorted_working_model_list:
      working_model_list.append(m)

    # Go through all the working models and cross them with other models to
    #  optimize...Then take all the best and cross...

    best_score,best_model=sorted_working_model_list[0]
    found=True
    cycle=0
    while found:
      cycle+=1
      print("\nCYCLE %d current best is %7.3f\n" %(
        cycle,best_model.get_score()), file=out)
      found=False
      sorted_working_model_list=[]
      new_best=best_model
      id=0
      for working_model in working_model_list:
        id+=1
        others=[]
        for m in working_model_list:
          if not working_model==m:  others.append(m)
        new_working_model=working_model.optimize_with_others(others=others)
        if not new_working_model:
          print()
          continue
        aa=[new_working_model.get_score(),new_working_model]
        if not aa in sorted_working_model_list:
          sorted_working_model_list.append(aa)
      if not sorted_working_model_list:
         break # nothing to do

      sorted_working_model_list = sorted(sorted_working_model_list,
         key = lambda wm: wm[0], reverse = True)
      sorted_working_model_list=sorted_working_model_list[:max_keep]

      new_working_score,new_working_model=sorted_working_model_list[0]
      if new_working_score>best_model.get_score():
        best_model=new_working_model
        found=True
        if verbose:
          print("NEW BEST SCORE: %7.2f" %(best_model.get_score()), file=out)
          best_model.show_summary(out=out)

    print("\nDONE... best is %7.3f\n" %(
        best_model.get_score()), file=out)

    # Create composite of this chain

    # Note residue values. We are going to pick each residue from one of
    # the models
    for model in ph.models():
      for chain in model.chains():
        if chain.id != chain_id: continue
        residue_list=[]
        for rg in chain.residue_groups():
          residue_list.append(rg.resseq)
    residue_list.sort()
    assert len(best_model.source_list)==len(residue_list)
    from mmtbx.secondary_structure.find_ss_from_ca import remove_ter_or_break
    for i in range(len(residue_list)):
      atom_selection=get_atom_selection(model_id=best_model.source_list[i],
        resseq_sel=residue_list[i])
      asc=ph.atom_selection_cache()
      sel=asc.selection(string = atom_selection)
      sel_hierarchy=ph.select(sel)
      sel_hierarchy = remove_ter_or_break(sel_hierarchy)
      sel_ph_list.append(sel_hierarchy)
  from iotbx.pdb.utils import add_hierarchies
  pdb_hierarchy = remove_ter_or_break(add_hierarchies(sel_ph_list,
    create_new_chain_ids_if_necessary = False))

  if pdb_out:
    pdb_out = pdb_hierarchy.write_pdb_or_mmcif_file(target_filename = pdb_out,
      crystal_symmetry = crystal_symmetry)
    print("Final model is in: %s\n" %(pdb_out))

  return pdb_hierarchy, pdb_out

if   (__name__ == "__main__"):
  args=sys.argv[1:]
  # Get the parameters
  params=get_params(args=args)
  run(params=params)
