from __future__ import division
import iotbx.phil
import iotbx.ccp4_map
from cctbx import crystal
from cctbx import maptbx
from libtbx.utils import Sorry
import sys, os
from cctbx.array_family import flex

master_phil = iotbx.phil.parse("""

  input_files {

    ccp4_map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    ncs_file = None
      .type = path
      .help = File with NCS information (typically point-group NCS with \
               the center specified). Typically in  PDB format. \
              Can also be a .ncs_spec file from phenix.
      .short_caption = NCS info file

    pdb_file = None
      .type = path
      .help = Optional PDB file matching ccp4_map_file to be offset

  }

  output_files {

    shifted_map_file = shifted_map.ccp4
      .type = path
      .help = Input map file shifted to new origin. Only written if a shift is\
                applied.
      .short_caption = Shifted map file

    shifted_pdb_file = shifted_pdb.pdb
      .type = path
      .help = Input pdb file shifted to new origin. Only written if a shift is\
                applied.
      .short_caption = Shifted pdb file

    shifted_ncs_file = shifted_ncs.ncs_spec 
      .type = path
      .help = NCS information shifted to new origin.  Only written if a shift \
         is applied.
      .short_caption = Output NCS info file


    output_map_directory =  None
      .type = path
      .help = Directory where output maps are to be written \
                applied.
      .short_caption = Output map directory


    box_map_file = box_map_au.ccp4
      .type = path
      .help = Output map file with one NCS asymmetric unit, cut out box
      .short_caption = Box NCS map file

    box_mask_file = box_mask_au.ccp4
      .type = path
      .help = Output mask file with one NCS asymmetric unit, cut out box
      .short_caption = Box NCS mask file

    box_buffer = 5
      .type = int
      .help = Buffer (grid units) around NCS asymmetric unit in box_mask and map
      .short_caption = Box buffer size

    au_map_output_file = shifted_map_au.ccp4
      .type = path
      .help = Output map file with one NCS asymmetric unit included
      .short_caption = Output NCS map file

    au_mask_output_file = shifted_mask_au.ccp4
      .type = path
      .help = Output mask file with one NCS asymmetric unit marked
      .short_caption = Output NCS mask file

    write_intermediate_maps = False
      .type = bool
      .help = Write out intermediate maps and masks for visualization
      .short_caption = Write intermediate maps

    write_output_maps = True
      .type = bool
      .help = Write out maps
      .short_caption = Write maps

    remainder_map_file = remainder_map.ccp4
      .type = path
      .help = output map file with remainder after initial regions identified
      .short_caption = Output remainder map file

    output_info_file = None
      .type = path
      .help = Output pickle file with information about map and masks
      .short_caption = Output pickle file
  }

  crystal_info {
     space_group = None
       .type = space_group
       .short_caption = Space Group
       .help = Space group (normally read from the data file)
       .style = bold
     unit_cell = None
       .type = unit_cell
       .short_caption = Unit Cell
       .help = Unit Cell (normally read from the data file)
       .style = bold
     seq_file = None
       .type = path
       .short_caption = Sequence file
       .help = Sequence file (unique chains only,  \
               1-letter code, chains separated by \
               blank line or greater-than sign.)  \
               Can have chains that are DNA/RNA/protein and\
               all can be present in one file.
     chain_type = *None PROTEIN RNA DNA
       .type = choice
       .short_caption = Chain type
       .help = Chain type. Determined automatically from sequence file if \
               not given.
  }

  segmentation {
    density_threshold = None
      .type = float
      .short_caption = Density threshold
      .help = Threshold density for identifying regions of density. \
             Applied after normalizing the density in the region of \
             the molecule to an rms of 1 and mean of zero.
    remove_bad_regions = True
      .type = bool
      .short_caption = Remove bad regions
      .help = Remove regions that are bigger than the NCS asymmetric unit.\

    split_if_possible = True
      .type = bool
      .short_caption = Split regions if mixed
      .help = Split regions that are split in some NCS copies.\
              If None, split if most copies are split.
    write_all_regions = False
      .type = bool
      .short_caption = Write all regions
      .help = Write all regions to ccp4 map files.

    max_per_au = None
      .type = int
      .short_caption = Max regions in au
      .help = Maximum number of regions to be kept in the NCS asymmetric unit

    min_ratio = 0.1
      .type = float
      .short_caption = Minimum ratio to keep
      .help = Minimum ratio of region size to maximum to keep it

    max_ratio_to_target = 3
      .type = float
      .help = Maximum ratio of grid points in top region to target
      .short_caption = Max ratio to target

    min_volume = 10
      .type = int
      .help = Minimum region size to consider (in grid points)
      .short_caption = Minimum region size

    max_overlap_fraction = 0.01
      .type = float
      .short_caption = Max overlap
      .help = Maximum fractional overlap allowed to density in another \
              asymmetric unit.

    iterate_with_remainder = True
      .type = bool
      .short_caption = Iterate
      .help = Iterate looking for regions based on remainder from first analysis

    expand_size = None
      .type = int
      .help = Grid points to expand size of regions when excluding for next \
               round. If None, set to approx number of grid points to get \
               expand_target below
      .short_caption = Expand size

    expand_target = 1.5
      .type = float
      .help = Target expansion of regions (A)
      .short_caption = Expand target


  }
   control {
      verbose = False
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output
      njump = 2
        .type = int
        .help = Take every njump'th grid point in each direction in analyses
        .short_caption = Grid skip interval
   }
""", process_includes=True)
master_params = master_phil

class info_object:
  def __init__(self,
      ncs_obj=None,
      ncs_group_list=None,
      origin_shift=None,
      edited_volume_list=None,
      region_range_dict=None,
      selected_regions=[],
      ncs_related_regions=[],
      map_files_written=[],
      bad_region_list=None,
      region_centroid_dict=None,
      original_id_from_id=None,
      remainder_id_dict=None,  # dict relating regions in a remainder object to
    ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

class ncs_group_object:
  def __init__(self,
      ncs_obj=None,
      ncs_group_list=None,
      edited_mask=None,
      origin_shift=None,
      edited_volume_list=None,
      region_range_dict=None,
      selected_regions=[],
      ncs_related_regions=[],
      map_files_written=[],
      bad_region_list=None,
      region_centroid_dict=None,
      region_scattered_points_dict=None,
      co=None,
      conn_obj=None,
      original_id_from_id=None,
      remainder_id_dict=None,  # dict relating regions in a remainder object to
                               #  those in the original map
         ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

  def as_info_object(self):
    return info_object(
      ncs_obj=self.ncs_obj,
      ncs_group_list=self.ncs_group_list,
      origin_shift=self.origin_shift,
      edited_volume_list=self.edited_volume_list,
      region_range_dict=self.region_range_dict,
      selected_regions=self.selected_regions,
      ncs_related_regions=self.ncs_related_regions,
      bad_region_list=self.bad_region_list,
      region_centroid_dict=self.region_centroid_dict,
      original_id_from_id=self.original_id_from_id,
      map_files_written=self.map_files_written,
     )
  def set_selected_regions(self,selected_regions):
    from copy import deepcopy
    self.selected_regions=deepcopy(selected_regions)

  def set_ncs_related_regions(self,ncs_related_regions):
    from copy import deepcopy
    self.ncs_related_regions=deepcopy(ncs_related_regions)

  def set_map_files_written(self,map_files_written):
    from copy import deepcopy
    self.map_files_written=deepcopy(map_files_written)

def write_ccp4_map(crystal_symmetry, file_name, map_data):
  iotbx.ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=crystal_symmetry.unit_cell(),
      space_group=crystal_symmetry.space_group(),
      map_data=map_data.as_double(),
      labels=flex.std_string([""]))

def set_up_xrs(crystal_symmetry=None):  # dummy xrs to write out atoms

  lines=["ATOM     92  SG  CYS A  10       8.470  28.863  18.423  1.00 22.05           S"] # just a random line to set up x-ray structure
  import iotbx.pdb
  from cctbx.array_family import flex
  from cctbx import xray
  pdb_inp=iotbx.pdb.input(source_info="",lines=lines)
  xrs = pdb_inp.xray_structure_simple(crystal_symmetry=crystal_symmetry)
  scatterers = flex.xray_scatterer()
  return xrs,scatterers

  """from cctbx import xray
  scatterers.append( xray.scatterer(scattering_type="O", label="O",
    site=xyz_frac, u=0.1, occupancy=1.0))
  """

def write_xrs(xrs=None,scatterers=None,file_name="atoms.pdb"):
  from cctbx import xray
  xrs = xray.structure(xrs, scatterers=scatterers)
  text=xrs.as_pdb_file()
  f=open(file_name,'w')
  print >>f,text
  f.close()
  print "Atoms written to %s" %file_name

def get_params(args,out=sys.stdout):
  import mmtbx.utils
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_phil)
  params = inputs.params.extract()
  print >>out,"\nSegment_and_split_map\n"
  print >>out,"Command used: %s\n" %(
   " ".join(['segment_and_split_map']+args))
  master_params.format(python_object=params).show(out=out)

  # PDB file
  if params.input_files.pdb_file and not inputs.pdb_file_names:
    inputs.pdb_file_names=[params.input_files.pdb_file]
  if inputs.pdb_file_names:
    params.input_files.pdb_file=inputs.pdb_file_names[0]
  print >>out,"\nInput PDB file: %s\n" %(params.input_files.pdb_file)
  if params.input_files.pdb_file:
    pdb_inp = iotbx.pdb.input(file_name=params.input_files.pdb_file)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
  else:
    pdb_hierarchy=None

  if inputs.ccp4_map:
    params.input_files.ccp4_map_file=inputs.ccp4_map_file_name
  else:
    if params.input_files.ccp4_map_file:
      from iotbx import ccp4_map
      inputs.ccp4_map=iotbx.ccp4_map.map_reader(
       file_name=params.input_files.ccp4_map_file)
  if not inputs.ccp4_map:
    raise Sorry("Need ccp4 map")

  map_data=inputs.ccp4_map.map_data()
  crystal_symmetry=crystal.symmetry(inputs.ccp4_map.unit_cell().parameters(),
    inputs.ccp4_map.space_group_number)

  shift_needed = not \
      (map_data.focus_size_1d() > 0 and map_data.nd() == 3 and
       map_data.is_0_based())

  a,b,c = crystal_symmetry.unit_cell().parameters()[:3]
  N_ = map_data.all()
  O_ =map_data.origin()
  sx,sy,sz = (a/N_[0])*O_[0], (b/N_[1])*O_[1], (c/N_[2])*O_[2]
  print >>out,"Origin for map is at (%8.2f,%8.2f,%8.2f)" % (-sx,-sy,-sz)
  print >>out,"Cell dimensions of this map are: (%8.2f,%8.2f,%8.2f)" % (a,b,c)
  if shift_needed:
    if(not crystal_symmetry.space_group().type().number() in [0,1]):
        raise RuntimeError("Not implemented")
    origin_shift=[-sx,-sy,-sz]
    print >>out, "Moving origin to (0,0,0)"
    print >>out,"Adding (%8.2f,%8.2f,%8.2f) to all coordinates\n"%(sx,sy,sz)

    map_data=map_data.shift_origin()
  else:
    origin_shift=None

  if params.segmentation.expand_size is None:

    abc = crystal_symmetry.unit_cell().parameters()[:3]
    N_ = map_data.all()
    nn=0.
    for i in xrange(3):
      delta=abc[i]/N_[i]
      nn+=params.segmentation.expand_target/delta
    nn=max(1,int(0.5+nn/3.))
    params.segmentation.expand_size=nn
    print >>out,"Expand size (grid units): %d (about %4.1f A) " %(
      nn,nn*abc[0]/N_[0])

  return params,crystal_symmetry,map_data,origin_shift,pdb_hierarchy

def get_ncs(params,out=sys.stdout):
  file_name=params.input_files.ncs_file
  if not file_name: # No ncs supplied...use just 1 ncs copy..
    from mmtbx.ncs.ncs import ncs
    ncs_object=ncs()
    ncs_object.set_unit_ncs()
    #ncs_object.display_all(log=out)
  elif not os.path.isfile(file_name):
    raise Sorry("The ncs file %s is missing" %(file_name))
  else: # get the ncs
    from mmtbx.ncs.ncs import ncs
    ncs_object=ncs()
    try: # see if we can read biomtr records
      pdb_inp=iotbx.pdb.input(file_name=file_name)
      ncs_object.ncs_from_pdb_input_BIOMT(pdb_inp=pdb_inp,log=out)
    except Exception,e: # try as regular ncs object
      ncs_object.read_ncs(file_name=file_name,log=out)
    #ncs_object.display_all(log=out)
    print >>out,"\nTotal of %d NCS operators read\n" %(
      ncs_object.max_operators())
  if not ncs_object or ncs_object.max_operators()<1:
    raise Sorry("Need ncs information from an ncs_info file")
  return ncs_object

def score_threshold(threshold=None,
     sorted_by_volume=None,n_residues=None,
     ncs_copies=None,
     target_fraction=None,
     solvent_fraction=None,
     map_data=None,
     residues_per_region=50.,
     min_volume=None,
     min_ratio=None,
     max_ratio_to_target=None,
     weight_score_grid_points=1.,
     weight_score_ratio=1.0,
     weight_near_one=0.1,
     minimum_ratio_of_ncs_copy_to_first=0.5,
     out=sys.stdout):

   # We want about 1 region per 50-100 residues for the biggest region.
   # One possibility is to try to maximize the median size of the N top
   # regions, where N=number of expected regions= n_residues/residues_per_region

   # Also note we have an idea how big a region should be (how many
   # grid points) if we make an assumption about the fractional volume that
   # should be inside a region compared to the total volume of protein/nucleic
   # acid in the region...this gives us target_in_top_regions points.
   # So using this, make the median size as close to target_in_top_regions as
   # we can.
   grid_points=map_data.size()
   expected_regions=max(1,int(0.5+n_residues/residues_per_region))

   target_in_all_regions=float(grid_points)*target_fraction*(1-solvent_fraction)
   target_in_top_regions=target_in_all_regions/expected_regions

   nn=len(sorted_by_volume)-1 # first one is total
   ok=True
   if nn < ncs_copies:
     ok=False #return  # not enough


   v1,i1=sorted_by_volume[1]
   if v1 < min_volume:
     ok=False #return

   if v1 > max_ratio_to_target*target_in_top_regions:
     ok=False #return

   # there should be about ncs_copies copies of each size region if ncs_copies>1
   if ncs_copies>1:
     v2,i2=sorted_by_volume[max(1,min(ncs_copies,nn))]
     score_ratio=v2/v1  # want it to be about 1
     if score_ratio < minimum_ratio_of_ncs_copy_to_first:
       ok=False #return  # not allowed
   else:
     score_ratio=1.0 # for ncs_copies=1

   nn2=min(nn,max(1,(len(sorted_by_volume[1:expected_regions])+1)//2))
   median_number,iavg=sorted_by_volume[nn2]

   # number in each region should be about target_in_top_regions
   if median_number > target_in_top_regions:
     score_grid_points=target_in_top_regions/max(1.,median_number)
   else:
     score_grid_points=median_number/target_in_top_regions

   if threshold>1.:
     score_near_one=1./threshold
   else:
     score_near_one=threshold

   overall_score=(
     (weight_score_ratio*score_ratio+
     weight_score_grid_points*score_grid_points+
     weight_near_one*score_near_one
       ) /
     (weight_score_ratio+weight_score_grid_points+weight_near_one))

   half_expected_regions=max(1,(1+expected_regions)//2)
   if ok and v1 >= target_in_top_regions/2 and \
        len(sorted_by_volume)>half_expected_regions:
     ratio=sorted_by_volume[half_expected_regions][0]/v1
     last_volume=sorted_by_volume[half_expected_regions][0]
     if ratio >=min_ratio and \
         last_volume>=min_volume:
       has_sufficient_regions=True
     else:
       has_sufficient_regions=False
   else:
       has_sufficient_regions=False
       ratio=0.

   print >>out,\
    " %5.2f   %5d     %4d    %5d     %5d     %6.3f   %5s    %5.3f" %(
       threshold,target_in_top_regions,expected_regions,
       v1,median_number,ratio,has_sufficient_regions,overall_score)

   if ok:
     return overall_score,has_sufficient_regions
   else:
     return


def choose_threshold(params,map_data=None,
     target_fraction=0.2,
     solvent_fraction=None,
     n_residues=None,
     ncs_copies=None,
     out=sys.stdout):

  best_threshold=None
  best_score=None

  print >>out,"\nChecking possible cutoffs for region identification"
  used_ranges=[]
  print >>out,\
    "Threshold  Target    N     Biggest   Median     Ratio    OK    Score"

  # Assume any threshold that is lower than a threshold that gave a non-zero value
  #  and is zero is an upper bound on the best value.  Same the other way around
  upper_bound=100
  lower_bound=0.01
  best_nn=None
  for n_range_low,n_range_high in [[-16,4],[-32,16],[-64,80]]:
    last_score=None
    for nn in xrange(n_range_low,n_range_high+1):
      if nn in used_ranges: continue
      used_ranges.append(nn)
      threshold=0.95**nn
      if threshold < lower_bound or threshold > upper_bound: 
        continue
      co = maptbx.connectivity(map_data=map_data.deep_copy(), 
         threshold=threshold)
      z = zip(co.regions(),range(0,co.regions().size()))
      sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
      if len(sorted_by_volume)<2:
        info=None
      else: 
        info=score_threshold(
         threshold=threshold,
         sorted_by_volume=sorted_by_volume,
         target_fraction=target_fraction,
         solvent_fraction=solvent_fraction,
         min_volume=params.segmentation.min_volume,
         min_ratio=params.segmentation.min_ratio,
         max_ratio_to_target=params.segmentation.max_ratio_to_target,
         ncs_copies=ncs_copies,
         n_residues=n_residues,
         map_data=map_data,
         out=out)
      if not info or not info[0]: # zero value
        if best_threshold and threshold >best_threshold: # new upper bound
           upper_bound=threshold
        elif best_threshold and threshold <best_threshold: # new upper bound
           lower_bound=threshold
        continue

      score,has_sufficient_regions=info
      if score and ( best_score is None or score > best_score):
        best_threshold=threshold
        best_score=score

  if params.segmentation.density_threshold is not None: # use it
     best_threshold=params.segmentation.density_threshold
     print >>out,"\nUsing input threshold of %5.2f " %(best_threshold)
  if best_threshold is None:
    raise Sorry("Threshold not successfully identified")

  print >>out,"\nBest threshold: %5.2f\n" %(best_threshold)
  return best_threshold

def get_connectivity(params,
     map_data=None,
     ncs_object=None,
     solvent_fraction=None,
     target_fraction=0.2,
     n_residues=None,
     ncs_copies=None,
     out=sys.stdout):
  print >>out,"\nGetting connectivity"

  # get map data and normalize to SD of the part that is not solvent
  sd=map_data.sample_standard_deviation()
  scaled_sd=sd/(1-solvent_fraction)**0.5
  map_data=(map_data-map_data.as_1d().min_max_mean().mean)/scaled_sd

  # Try connectivity at various thresholds
  # Choose one that has about the right number of grid points in top regions
  threshold=choose_threshold(params,
     map_data=map_data,
     n_residues=n_residues,
     ncs_copies=ncs_copies,
     target_fraction=target_fraction,
     solvent_fraction=solvent_fraction,
     out=out)

  co = maptbx.connectivity(map_data=map_data, threshold=threshold)
  conn = co.result()
  z = zip(co.regions(),range(0,co.regions().size()))
  sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)

  min_b, max_b = co.get_blobs_boundaries_tuples() # As grid points, not A

  cntr=0
  v1=sorted_by_volume[1][0]
  n_use=0
  for p in sorted_by_volume[1:]:
    cntr+=1
    v,i=p
    if(v<params.segmentation.min_volume): break
    if(v/v1 <params.segmentation.min_ratio): break
    """
    print >>out,\
    "Region %3d (%3d)  volume:%5d  X:%6d - %6d   Y:%6d - %6d  Z:%6d - %6d "%(
     cntr,i,v,
     min_b[i][0],max_b[i][0],
     min_b[i][1],max_b[i][1],
     min_b[i][2],max_b[i][2])
    """
    n_use=cntr
  return co,conn,sorted_by_volume,min_b,max_b

def get_volume_of_seq(text,vol=None,chain_type=None,out=sys.stdout):
  chain_type,n_residues=guess_chain_type(text,chain_type=chain_type,out=out)
  if chain_type is None and n_residues is None:
    return None,None

  if chain_type=='PROTEIN':
    mw_residue=110.0  # from $CDOC/matthews.doc
    density_factor=1.23   # 1.66/DENSITY-OF-PROTEIN=1.66/1.35
  else:
    mw_residue=330.0  # guess for DNA
    density_factor=1.15   # 1.66/DENSITY-OF-DNA=1.66/1.45
  return len(text)*density_factor*mw_residue,len(text)

def create_rna_dna(cns_dna_rna_residue_names):
  dd={}
  for key in cns_dna_rna_residue_names.keys():
    dd[cns_dna_rna_residue_names[key]]=key
  return dd

def guess_chain_type(text,chain_type=None,out=sys.stdout):
    if chain_type in ['None',None]:
      test_types=['RNA','PROTEIN']
    else:
      test_types=[chain_type]
    from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter
    from iotbx.pdb import cns_dna_rna_residue_names
    three_letter_rna_dna_given_one_letter=create_rna_dna(
       cns_dna_rna_residue_names)

    residues_dict={}
    most_residues=0
    best="RNA" # start with DNA so we separate RNA/DNA first, then protein
    found_non_gly_non_ala=False
    text=text.upper()
    for test_chain_type in test_types:
      residues_dict[test_chain_type]=0
      seq_text=""
      for c in text:
        if test_chain_type in ['RNA','DNA']:
          three_char=three_letter_rna_dna_given_one_letter.get(c,"")
        else:
          three_char=three_letter_given_one_letter.get(c,"")
        if not three_char:continue
        seq_text+=c
        if test_chain_type=='PROTEIN' and three_char not in ['GLY','ALA']:
          found_non_gly_non_ala=True

      residues_dict[test_chain_type]=len(seq_text)
      if residues_dict[test_chain_type]>most_residues:
        most_residues=residues_dict[test_chain_type]
        best=test_chain_type

    if best in ['RNA','DNA'] and \
        residues_dict.get('PROTEIN',0)==most_residues and \
        not found_non_gly_non_ala:  # guess protein here
      best='PROTEIN'
    if most_residues<1:
      return None,None
    print >>out,"\nChain type for chain with %d residues set to %s" %(
       most_residues,best)
    return best,most_residues

def get_solvent_fraction(params,crystal_symmetry=None,
     ncs_object=None,out=sys.stdout):
  map_volume=crystal_symmetry.unit_cell().volume()
  ncs_copies=ncs_object.max_operators()
  if not params.crystal_info.seq_file:
    raise Sorry("Please specify a sequence file with seq_file=myseq.seq")
  elif not os.path.isfile(params.crystal_info.seq_file):
    raise Sorry(
     "The sequence file '%s' is missing." %(params.crystal_info.seq_file))
  seq_as_string=open(params.crystal_info.seq_file).read()
  seq_as_string=">\n"+seq_as_string  # so it always starts with >
  seq_as_string=seq_as_string.replace("\n\n","\n>\n") # blank lines are like >
  spl=seq_as_string.split(">")
  volume_of_chains=0.
  n_residues=0
  for s in spl:
    if not s: continue
    ss="".join(s.splitlines()[1:])
    volume,nres=get_volume_of_seq(ss,vol=map_volume,
      chain_type=params.crystal_info.chain_type,out=out)
    if volume is None: continue
    volume_of_chains+=volume
    n_residues+=nres
  volume_of_molecules=volume_of_chains*ncs_copies
  n_residues=n_residues*ncs_copies
  solvent_fraction=1.-(volume_of_molecules/map_volume)
  solvent_fraction=max(0.001,solvent_fraction)
  print >>out, \
    "Cell volume: %.1f  NCS copies: %d   Volume of unique chains: %.1f" %(
     map_volume,ncs_copies,volume_of_chains)
  print >>out,\
    "Total residues: %d  Volume of all chains: %.1f  Solvent fraction: %.3f "%(
       n_residues,volume_of_molecules,solvent_fraction)
  return solvent_fraction,n_residues

def mask_for_region(id=None,conn_obj=None):
   region_mask_i = conn_obj.deep_copy()
   s = (region_mask_i==id)
   region_mask_i = region_mask_i.set_selected(s,1)
   region_mask_i = region_mask_i.set_selected(~s,0)
   return region_mask_i

def region_is_in_multiple_ncs_au(id=None,
        conn_obj=None,ncs_obj=None,
        max_overlap=None,out=sys.stdout):
   # figure out if any points in this region are in more than one ncs au
   # map all points in the region to ncs copies and see if they are then in
   #  the same region
   # value_at_closest_grid_point
   region_mask_i=mask_for_region(id=id,conn_obj=conn_obj)
   return False

def top_key(dd):
  if not dd:
    return None,None
  elif len(dd.keys())==1:
    return dd.keys()[0],dd[dd.keys()[0]]
  else:
    best_key=None
    best_n=None
    for key in dd.keys():
      if not best_n or dd[key] > best_n:
        best_n=dd[key]
        best_key=key
    return best_key,best_n

def choose_max_regions_to_consider(params,
    sorted_by_volume=None,
    ncs_copies=None):

  max_per_au=params.segmentation.max_per_au
  min_ratio=params.segmentation.min_ratio
  min_volume=params.segmentation.min_volume
  # sort and eliminate regions with few points and those at end of list
  if len(sorted_by_volume)<2:
    return 0
  max_grid_points=sorted_by_volume[1][0]
  cntr=0
  for p in sorted_by_volume[1:]:
    cntr+=1
    if max_per_au and (cntr>max_per_au*ncs_copies):
      cntr-=1
      break
    v,i=p  # v=volume in grid points, i=id
    if v/max_grid_points<min_ratio or v < min_volume:
      cntr-=1
      break
  return cntr

def get_edited_mask(sorted_by_volume=None,
    conn_obj=None,max_regions_to_consider=None,
    out=sys.stdout):
  origin=list(conn_obj.accessor().origin())
  all=list(conn_obj.accessor().all())
  conn_obj.accessor().show_summary(out)
  edited_mask=conn_obj.deep_copy()
  first=True
  edited_volume_list=[]
  original_id_from_id={}
  for i in xrange(1,max_regions_to_consider+1):
    v,id=sorted_by_volume[i]
    original_id_from_id[i]=id
    edited_volume_list.append(v)
    s = (conn_obj==id)
    if first:
      edited_mask=edited_mask.set_selected(~s,0)
      first=False
    edited_mask=edited_mask.set_selected(s,i)   # edited mask has ID of
         # regions, labeled in decreasing size from 1 to max_regions_to_consider
  return edited_mask,edited_volume_list,original_id_from_id

def accumulate(region_range_dict,region_centroid_dict,region_n_dict,
          region_scattered_points_dict,
          xyz_cart=None,id=None):
  if not id in region_n_dict.keys():
    region_n_dict[id]=0
    region_scattered_points_dict[id]=flex.vec3_double()
    region_centroid_dict[id]=list(xyz_cart)
    region_range_dict[id]=[ [None,None],[None,None],[None,None] ]
  region_n_dict[id]+=1
  region_scattered_points_dict[id].append(xyz_cart)
  a=list(xyz_cart)
  region_centroid_dict[id][0]+=a[0]
  region_centroid_dict[id][1]+=a[1]
  region_centroid_dict[id][2]+=a[2]
  for pair,x in zip(region_range_dict[id],xyz_cart):
    if pair[0] is None or pair[0]>x: pair[0]=x
    if pair[1] is None or pair[1]<x: pair[1]=x

def choose_subset(a,target_number=1):
  new_array=flex.vec3_double()
  assert type(new_array)==type(a)
  n=a.size()
  njump=max(1,n//target_number)
  i=0
  for x in a:
    if i%njump==0 or i==n-1:
     new_array.append(x)
    i+=1
  return new_array

def analyze(region_range_dict,region_centroid_dict,region_n_dict,
     region_scattered_points_dict,
     out=sys.stdout):
  #print >>out,"Centroid and ranges for regions:"
  for id in region_n_dict.keys():
    for i in xrange(3):
      region_centroid_dict[id][i]=region_centroid_dict[id][i]/region_n_dict[id]
    rc=region_centroid_dict[id]
    region_scattered_points_dict[id]=choose_subset(
      region_scattered_points_dict[id],
      target_number=min(30,max(region_n_dict[id]/10,10)))

def run_get_duplicates_and_ncs(
   ncs_obj=None,
   conn_obj=None,
   crystal_symmetry=None,
   edited_mask=None,
   njump=None,
   unit_cell=None,
   max_regions_to_consider=None,
   out=sys.stdout,
   ):
  for njump in xrange(njump+1,0,-1):

    duplicate_dict,equiv_dict,equiv_dict_ncs_copy_dict,\
        region_range_dict,region_centroid_dict,region_scattered_points_dict=\
      get_duplicates_and_ncs(
        ncs_obj=ncs_obj,
        conn_obj=conn_obj,
        edited_mask=edited_mask,
        unit_cell=unit_cell,
        max_regions_to_consider=max_regions_to_consider,
        njump=njump,
        out=out)

    # check that we have region_centroid for all values
    complete=True
    missing=[]
    for i in xrange(1,max_regions_to_consider+1):
      if not i in region_centroid_dict.keys():
         complete=False
         missing.append(i)
    if complete:
       return duplicate_dict,equiv_dict,equiv_dict_ncs_copy_dict,\
        region_range_dict,region_centroid_dict,njump,region_scattered_points_dict
  raise Sorry("Cannot find region-centroid for all regions? Missing: %s" %(
      missing))


def get_duplicates_and_ncs(
   ncs_obj=None,
   conn_obj=None,
   edited_mask=None,
   njump=None,
   crystal_symmetry=None,
   unit_cell=None,
   max_regions_to_consider=None,
   write_atoms=None, # ID of region to write dummy atoms for (if any)
   out=sys.stdout,
   ):

  ncs_group=ncs_obj.ncs_groups()[0]
  duplicate_dict={}  # keyed by id, number of duplicates for that region
  equiv_dict={}  # equiv_dict[id][other_id]=number_of points other_id matches
                 #  id through an ncs relationship
  equiv_dict_ncs_copy_dict={}

  region_range_dict={} # keyed by region in edited_mask; range for x, y, z
  region_centroid_dict={} # median values of x,y,z
  region_n_dict={}  # count of points used by region (differs from volume due
     # to the njump parameter. Should be about volume/njump**3)
  region_scattered_points_dict={} # some points in each region

  cntr=0
  from scitbx.math import matrix

  origin=list(conn_obj.accessor().origin())
  all=list(conn_obj.accessor().all())

  if write_atoms:
    xrs,scatterers=set_up_xrs(crystal_symmetry=crystal_symmetry)

  for i in xrange (origin[0],all[0],njump):
    for j in xrange (origin[1],all[1],njump):
      for k in xrange (origin[2],all[2],njump):
        id=edited_mask[i,j,k]
        if id <1: continue

        # What regions are ncs-related points in?
        xyz_frac=(i/all[0],j/all[1],k/all[2])
        xyz_cart=unit_cell.orthogonalize(xyz_frac)
        accumulate(region_range_dict,region_centroid_dict,region_n_dict,
          region_scattered_points_dict,
          xyz_cart=xyz_cart,id=id)
        if write_atoms and id==write_atoms:
          from cctbx import xray
          scatterers.append( xray.scatterer(scattering_type="O", label="O",
           site=xyz_cart, u=0.1, occupancy=1.0))

        # Assume first one is identity
        n=0
        for t,r in zip(ncs_group.translations_orth()[1:],
            ncs_group.rota_matrices()[1:]):
          n+=1
          new_xyz_cart=r * matrix.col(xyz_cart) + t
          new_xyz_frac=unit_cell.fractionalize(new_xyz_cart)
          value=edited_mask.value_at_closest_grid_point(new_xyz_frac)
          if value==id:
            if not id in duplicate_dict.keys(): duplicate_dict[id]=0
            duplicate_dict[id]+=1
            break # only count once
          elif value>0:  # notice which one is matched
            if not id in equiv_dict.keys():
              equiv_dict[id]={}
              equiv_dict_ncs_copy_dict[id]={}
            if not value in equiv_dict[id].keys():
              equiv_dict[id][value]=0
              equiv_dict_ncs_copy_dict[id][value]={}
            equiv_dict[id][value]+=1
            if not n in equiv_dict_ncs_copy_dict[id][value].keys():
              equiv_dict_ncs_copy_dict[id][value][n]=0
            equiv_dict_ncs_copy_dict[id][value][n]+=1  # how many are ncs copy n
        cntr+=1

  if write_atoms:
    write_xrs(xrs=xrs,scatterers=scatterers,file_name="atoms.pdb")

  analyze(region_range_dict,region_centroid_dict,region_n_dict,
     region_scattered_points_dict,out=out)

  return duplicate_dict,equiv_dict,equiv_dict_ncs_copy_dict,\
      region_range_dict,region_centroid_dict,region_scattered_points_dict

def remove_bad_regions(params=None,
  duplicate_dict=None,
  edited_volume_list=None,
  out=sys.stdout):

  dups=duplicate_dict.keys()
  bad_region_list=[]
  if dups:
    print >>out,"\nRegions that span multiple NCS au:"
    dups.sort()
    for id in dups:
      print >>out,"ID: %d  Duplicate points: %d " %(
        id,duplicate_dict[id])

  new_sorted_by_volume=[]
  cntr=0
  region_list=[]
  region_volume_dict={}
  for i in xrange(len(edited_volume_list)):
    id=i+1
    cntr+=1
    v=edited_volume_list[i]
    if cntr in dups and \
       duplicate_dict[cntr]>params.segmentation.max_overlap_fraction*v \
        and params.segmentation.remove_bad_regions: #  mark it for not using
      bad_region_list.append(id)
    new_sorted_by_volume.append([v,id])
    region_list.append(id)
    region_volume_dict[id]=v
  if bad_region_list:
    print >>out,"Bad regions (excluded)",bad_region_list
  return region_list,region_volume_dict,new_sorted_by_volume,bad_region_list

def sort_by_ncs_overlap(matches,equiv_dict_ncs_copy_dict_id):
    sort_list=[]
    for id1 in matches:
      key,n=top_key(equiv_dict_ncs_copy_dict_id[id1]) # Take top ncs_copy
      sort_list.append([n,id1])
    sort_list.sort()
    sort_list.reverse()
    key_list=[]
    for n,id1 in sort_list:
      key_list.append(id1)
    return key_list


def get_ncs_equivalents(
    bad_region_list=None,
    region_list=None,
    region_volume_dict=None,
    equiv_dict=None,
    njump=None,
    ncs_copies=None,
    equiv_dict_ncs_copy_dict=None,
    min_coverage=.10,
    out=sys.stdout):

  equiv_dict_ncs_copy={}
  for id in region_list:
    if id in bad_region_list: continue
    match_dict=equiv_dict.get(id,{}) # which are matches
    matches=match_dict.keys()
    if not matches: continue
    key_list=sort_by_ncs_overlap(matches,equiv_dict_ncs_copy_dict[id])
    n_found=0
    for id1 in key_list:
      #     id matches id1 N=match_dict[id1]
      key,n=top_key(equiv_dict_ncs_copy_dict[id][id1]) # ncs_copy, n-overlap
      if n<min_coverage*region_volume_dict[id1]/(njump**3):
        break
      else:
        if not id in equiv_dict_ncs_copy.keys():equiv_dict_ncs_copy[id]={}
        equiv_dict_ncs_copy[id][id1]=key
        n_found+=1
        if n_found>=ncs_copies-1:
          break
 
  return equiv_dict_ncs_copy

  print >>out,"\nSets of NCS-related regions"
  keys=equiv_dict_ncs_copy.keys()
  keys.sort()
  used=[]
  for id in keys:
    #if id in used: continue
    others=equiv_dict_ncs_copy[id].keys()
    used+=others
    print >>out,"%d: " %(id),
    for id1 in others:
      key,n=top_key(equiv_dict_ncs_copy_dict[id][id1])
      print >>out,"%d:%d" %(id1,n),
    print >>out
  print >>out

def group_ncs_equivalents(
    region_list=None,
    region_volume_dict=None,
    equiv_dict_ncs_copy=None,
    ncs_copies=None,
    split_if_possible=None,
    out=sys.stdout):

  # equiv_dict_ncs_copy[id][id1]=ncs_copy
  # group together all the regions that are related to region 1...etc
  # if split_if_possible then skip all groups with multiple entries

  ncs_equiv_groups_as_list=[]
  ncs_equiv_groups_as_dict={}
  for id in region_list:
    equiv_group={}  #equiv_group[ncs_copy]=[id1,id2,id3...]
    equiv_group[0]=[id] # always
    for id1 in equiv_dict_ncs_copy.get(id,{}).keys():
      ncs_copy=equiv_dict_ncs_copy[id][id1]
      if not ncs_copy in equiv_group.keys(): equiv_group[ncs_copy]=[]
      equiv_group[ncs_copy].append(id1) # id1 is ncs_copy of id
    complete=True
    all_single=True
    equiv_group_as_list=[]
    total_grid_points=0
    for ncs_copy in xrange(ncs_copies): # goes 0 to ncs_copies-1
      equiv_group_as_list.append(equiv_group.get(ncs_copy,[]))
      if ncs_copy > 0 and \
          len(equiv_group.get(ncs_copy,[]))>1 and len(equiv_group.get(0,[]))==1:
        all_single=False
      for id in equiv_group.get(ncs_copy,[]):
        total_grid_points+=region_volume_dict[id]
      if not ncs_copy in equiv_group.keys():
        complete=False
    equiv_group_as_list.sort()
    if complete and \
        (not str(equiv_group_as_list) in ncs_equiv_groups_as_dict.keys() or
         total_grid_points>ncs_equiv_groups_as_dict[str(equiv_group_as_list)]) \
        and (all_single or (not split_if_possible)):
      ncs_equiv_groups_as_dict[str(equiv_group_as_list)]=total_grid_points
      ncs_equiv_groups_as_list.append([total_grid_points,equiv_group_as_list])

  ncs_equiv_groups_as_list.sort()
  ncs_equiv_groups_as_list.reverse()

  # Now remove any group that duplicates a previous group
  # 2015-11-07 allow a member to be in multiple groups though (for example
  #   one that spans several groups because it contains 2 region in other ncs
  #   copies)

  max_duplicates=ncs_copies-1  # don't let there be all duplicates
  ncs_group_list=[]
  used_list=[]
  print >>out,"All equiv groups:"
  used_regions=[]
  for total_grid_points,equiv_group_as_list in ncs_equiv_groups_as_list:
    duplicate=False
    for equiv_group in equiv_group_as_list:
      n_dup=0
      for x in equiv_group:
        if x in used_list:
          n_dup+=1
      if n_dup>max_duplicates:
        duplicate=True
    if not duplicate:
      print >>out,"NCS GROUP:",equiv_group_as_list,":",total_grid_points
      ncs_group_list.append(equiv_group_as_list)
      for equiv_group in equiv_group_as_list:
        for x in equiv_group: used_list.append(x)
  return ncs_group_list


def identify_ncs_regions(params,
     sorted_by_volume=None,
     co=None,
     conn_obj=None,
     ncs_obj=None,
     unit_cell=None,
     crystal_symmetry=None,
     origin_shift=None,
     out=sys.stdout):

  # 1.choose top regions to work with
  # 2.remove regions that are in more than one au of the NCS
  # 3.identify groups of regions that are related by NCS
  #  Also note the centers and bounds of each region

  # Choose number of top regions to consider

  max_regions_to_consider=choose_max_regions_to_consider(params,
    sorted_by_volume=sorted_by_volume,
    ncs_copies=ncs_obj.max_operators())

  print >>out,\
    "\nIdentifying NCS-related regions.Total regions to consider: %d" %(
    max_regions_to_consider)
  if max_regions_to_consider<1:
    print >>out,"\nUnable to identify any NCS regions"
    return None

  # Go through all grid points; discard if not in top regions
  #  Renumber regions in order of decreasing size

  edited_mask,edited_volume_list,original_id_from_id=get_edited_mask(
     sorted_by_volume=sorted_by_volume,
     conn_obj=conn_obj,
     max_regions_to_consider=max_regions_to_consider,out=out)
  # edited_mask contains re-numbered region id's

  # Identify duplicate and ncs relationships between regions
  # duplicate_dict[id]= number of duplicates for that region
  # equiv_dict[id][other_id]=number_of points other_id matches
                   #  id through an ncs relationship

  duplicate_dict,equiv_dict,equiv_dict_ncs_copy_dict,\
      region_range_dict,region_centroid_dict,params.control.njump,\
      region_scattered_points_dict=\
    run_get_duplicates_and_ncs(
      ncs_obj=ncs_obj,
      conn_obj=conn_obj,
      crystal_symmetry=crystal_symmetry,
      edited_mask=edited_mask,
      unit_cell=unit_cell,
      max_regions_to_consider=max_regions_to_consider,
      njump=params.control.njump,
      out=out)

  # Remove any bad regions
  region_list,region_volume_dict,new_sorted_by_volume,\
      bad_region_list=remove_bad_regions(
    params=params,
    duplicate_dict=duplicate_dict,
    edited_volume_list=edited_volume_list,
    out=out)

  # Identify groups of regions that are ncs-related
  # equiv_dict_ncs_copy[id][id1]=ncs_copy of id that corresponds to id1
  equiv_dict_ncs_copy=get_ncs_equivalents(
    region_list=region_list,
    bad_region_list=bad_region_list,
    region_volume_dict=region_volume_dict,
    equiv_dict=equiv_dict,
    ncs_copies=ncs_obj.max_operators(),
    njump=params.control.njump,
    equiv_dict_ncs_copy_dict=equiv_dict_ncs_copy_dict,
    out=out)

  # Group together regions that are ncs-related. Also if one ncs
  #   copy has 2 or more regions linked together, group the other ones.

  # each entry in ncs_group_list is a list of regions for each ncs_copy:
  #  e.g.,  [[8], [9, 23], [10, 25], [11, 27], [12, 24], [13, 22], [14, 26]]
  #  May contain elements that are in bad_region_list (to exclude later)
  ncs_group_list=group_ncs_equivalents(
    split_if_possible=params.segmentation.split_if_possible,
    ncs_copies=ncs_obj.max_operators(),
    region_volume_dict=region_volume_dict,
    region_list=region_list,
    equiv_dict_ncs_copy=equiv_dict_ncs_copy,
    out=out)

  ncs_group_obj=ncs_group_object(
     ncs_group_list=ncs_group_list,
     ncs_obj=ncs_obj,
     edited_mask=edited_mask,
     origin_shift=origin_shift,
     co=co,
     conn_obj=conn_obj,
     bad_region_list=bad_region_list,
     original_id_from_id=original_id_from_id,
     edited_volume_list=edited_volume_list,
     region_range_dict=region_range_dict,
     region_scattered_points_dict=region_scattered_points_dict,
     region_centroid_dict=region_centroid_dict)

  return ncs_group_obj

def get_center_list(regions,
    region_centroid_dict=None):
  center_list=[]
  for region in regions:
    center_list.append(region_centroid_dict[region])
  return center_list

def get_average_center(regions,
    region_centroid_dict=None):
  from copy import deepcopy
  center_list=get_center_list(regions,region_centroid_dict=region_centroid_dict)
  for region in regions:
    center_list.append(region_centroid_dict[region])
  average_center=deepcopy(center_list[0])
  if len(center_list)>1:
    for r in center_list[1:]:
      for i in xrange(3):
        average_center[i]+=r[i]
    for i in xrange(3):
      average_center[i]/=len(center_list)
  return average_center

def get_dist(r,s):
  dd=0.
  for i in xrange(3):
    dd+=(r[i]-s[i])**2
  return dd**0.5

def has_intersection(set1,set2):
  set1a=single_list(set1)
  set2a=single_list(set2)
  for x in set1a:
    if x in set2a:
      return True
  return False

def get_scattered_points_list(other_regions,
       region_scattered_points_dict=None):
  scattered_points_list=flex.vec3_double()
  for x in other_regions:
    scattered_points_list.extend(region_scattered_points_dict[x])
  return scattered_points_list


def get_closest_neighbor_rms(ncs_group_obj=None,selected_regions=None,
    target_scattered_points=None):
  # return rms closest distance of each region center to nearest other

  rms=0.
  rms_n=0.
  for x in selected_regions:
    test_centers=ncs_group_obj.region_scattered_points_dict[x]
    other_regions=remove_one_item(selected_regions,item_to_remove=x)
    other_center_list=get_scattered_points_list(other_regions,
       region_scattered_points_dict=ncs_group_obj.region_scattered_points_dict)
    if target_scattered_points:
      other_center_list.extend(target_scattered_points)
    dist=get_closest_dist(test_centers,other_center_list)
    if dist is not None:
      rms+=dist**2
      rms_n+=1.
  if rms_n>1:
    rms/=rms_n
  return rms**0.5


def get_rms(selected_regions=None,
    region_centroid_dict=None):
  # return rms distance of each region center from average of all others
  rms=0.
  rms_n=0.
  for x in selected_regions:
    other_regions=remove_one_item(selected_regions,item_to_remove=x)
    current_center=get_average_center(other_regions,
       region_centroid_dict=region_centroid_dict)
    test_center=region_centroid_dict[x]
    dist=get_dist(current_center,test_center)
    rms+=dist**2
    rms_n+=1.
  if rms_n>1:
    rms/=rms_n
  return rms**0.5

def single_list(list_of_lists):
  single=[]
  for x in list_of_lists:
    if type(x)==type([1,2,3]):
      single+=single_list(x)
    else:
      single.append(x)
  return single

def get_closest_dist(test_center,target_centers):
  closest_dist=None
  for target in target_centers:
    if type(test_center) in [type([1,2,3]),type(flex.vec3_double())]:
      dist=get_closest_dist(target,test_center)
    else:
      dist=get_dist(target,test_center)
    if closest_dist is None or dist<closest_dist:
      closest_dist=dist
  return closest_dist

def select_from_seed(starting_region,
      target_scattered_points=None,
      ncs_group_obj=None):
  from copy import deepcopy
  selected_regions=single_list(deepcopy(starting_region))
  # do not allow any region in ncs_group_obj.bad_region_list

  for ncs_group in ncs_group_obj.ncs_group_list: # try adding from each group
    best_ncs_set=None
    best_dist=None

    current_scattered_points_list=get_scattered_points_list(selected_regions,
       region_scattered_points_dict=ncs_group_obj.region_scattered_points_dict)
    if target_scattered_points:
      current_scattered_points_list.extend(target_scattered_points)
    for ncs_set in ncs_group: # pick the best ncs_set from this group
      if has_intersection(selected_regions,ncs_set) or \
         has_intersection(ncs_group_obj.bad_region_list,ncs_set):
        continue

      skip=False
      for x in selected_regions: # does any ncs copy of x overlap?
        if skip: continue
        for test_ncs_group in ncs_group_obj.ncs_group_list:
          if x in single_list(test_ncs_group):
            if has_intersection(test_ncs_group,ncs_set):
              skip=True
      if skip: continue

      # Get dist of this ncs_set (usually 1 region) to current center

      test_scattered_points=get_scattered_points_list(ncs_set,
       region_scattered_points_dict=ncs_group_obj.region_scattered_points_dict)

      dist=get_closest_dist(test_scattered_points,
        target_centers=current_scattered_points_list)
      if best_dist is None or dist<best_dist:
        best_dist=dist
        best_ncs_set=ncs_set
    if best_ncs_set is not None:
      selected_regions+=best_ncs_set

  rms=get_closest_neighbor_rms(ncs_group_obj=ncs_group_obj,
    target_scattered_points=target_scattered_points,
    selected_regions=selected_regions)

  return selected_regions,rms

def remove_one_item(input_list,item_to_remove=None):
  new_list=[]
  for item in input_list:
    if item != item_to_remove:
      new_list.append(item)
  return new_list

def get_ncs_related_regions(
    ncs_group_obj=None,
    selected_regions=None):

  ncs_related_regions=[]
  for ncs_group in ncs_group_obj.ncs_group_list:
    ids_in_group=single_list(ncs_group)
    for id in selected_regions:
      if id in ids_in_group:  # this group contains this selected id
        for i in ids_in_group:
          if (not i==id) and (not i in ncs_related_regions):
            ncs_related_regions.append(i)

  return ncs_related_regions

def all_elements_are_length_one(list_of_elements):
  for x in list_of_elements:
    if type(x)==type([1,2,3]):
      if len(x)!=1: return False
  return True

def select_regions_in_au(params,
     ncs_group_obj=None,
     target_scattered_points=None,
     out=sys.stdout):
  # Choose one region or set of regions from each ncs_group
  # Optimize closeness of centers...
  # If target scattered_points is supplied, include them as allowed target

  if not ncs_group_obj.ncs_group_list:
    return ncs_group_obj,[]

  if all_elements_are_length_one(ncs_group_obj.ncs_group_list):
    best_selected_regions=single_list(ncs_group_obj.ncs_group_list)
    best_rms=None
  else:

    best_selected_regions=None
    best_rms=None
    for starting_region in ncs_group_obj.ncs_group_list[0]: # try each au as seed
      if starting_region in ncs_group_obj.bad_region_list: continue # do not use
      selected_regions,rms=select_from_seed([starting_region],
        target_scattered_points=target_scattered_points,
        ncs_group_obj=ncs_group_obj)
      if not selected_regions:
        continue
      if best_rms is None or rms<best_rms:
        best_rms=rms
        best_selected_regions=selected_regions
        print >>out,"New best selected: rms: %7.1f: %s " %(
           rms,str(selected_regions))
    selected_regions=best_selected_regions
    if not selected_regions:
      raise Sorry("No NCS regions found from NCS groups")
    changing=True
    while changing:
      # Now see if replacing any regions with alternatives would improve it
      changing=False
      for x in selected_regions:
        starting_regions=remove_one_item(selected_regions,item_to_remove=x)
        new_selected_regions,rms=select_from_seed(starting_regions,
          target_scattered_points=target_scattered_points,
          ncs_group_obj=ncs_group_obj)
        if not new_selected_regions: continue
        if best_rms is None or rms<best_rms-0.00001:
          changing=True
          best_rms=rms
          best_selected_regions=new_selected_regions
          print >>out,"Optimized best selected: rms: %7.1f: %s " %(
            rms,str(selected_regions))

  selected_regions=best_selected_regions
  selected_regions.sort()

  print >>out,"\nFinal selected regions: ",
  for x in selected_regions:
    print >>out,x,
  print >>out

  # Identify scattered points for all selected regions:

  scattered_points=get_scattered_points_list(selected_regions,
     region_scattered_points_dict=ncs_group_obj.region_scattered_points_dict)

  # Identify ncs-related regions for all the selected regions
  ncs_related_regions=get_ncs_related_regions(
    ncs_group_obj=ncs_group_obj,
    selected_regions=selected_regions)
  print >>out,"NCS-related regions (not used): ",
  for x in ncs_related_regions:
     print >>out,x,
  print >>out
  ncs_group_obj.set_selected_regions(selected_regions)
  ncs_group_obj.set_ncs_related_regions(ncs_related_regions)

  return ncs_group_obj,scattered_points

def get_bool_mask_as_int(ncs_group_obj=None,mask_as_bool=None):
  mask_as_int=ncs_group_obj.edited_mask.deep_copy()
  s = (mask_as_bool==True)
  mask_as_int = mask_as_int.set_selected(s,1)
  mask_as_int = mask_as_int.set_selected(~s,0)
  return mask_as_int

def get_bool_mask_of_regions(ncs_group_obj=None,region_list=None,
    expand_size=None):
  s = (ncs_group_obj.edited_mask == -1)
  if region_list is None: region_list=[]
  for id in region_list:

    if not expand_size:
      s |= (ncs_group_obj.edited_mask==id)  # just take this region

    else:  # expand the size of the regions...use expand_mask which operates
         # on the original id numbers and uses the co
      bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand=ncs_group_obj.original_id_from_id[id],
        expand_size=expand_size)
      s |= (bool_region_mask== True)

  bool_mask = ncs_group_obj.co.expand_mask(id_to_expand=1,expand_size=1) # just to get bool mask
  bool_mask = bool_mask.set_selected(s,True)
  bool_mask = bool_mask.set_selected(~s,False)

  return bool_mask

def create_remaining_mask_and_map(params,
    ncs_group_obj=None,
    map_data=None,
    crystal_symmetry=None,
    out=sys.stdout):

  if not ncs_group_obj.selected_regions:
    print >>out,"No regions selected"
    return map_data


  # create new remaining_map containing everything except the part that
  # has been interpreted (and all points in interpreted NCS-related copies)

  bool_all_used=get_bool_mask_of_regions(ncs_group_obj=ncs_group_obj,
     region_list=ncs_group_obj.selected_regions+ncs_group_obj.ncs_related_regions,
     expand_size=params.segmentation.expand_size)
  map_data_remaining=map_data.deep_copy()
  s=(bool_all_used==True)
  map_data_remaining=map_data_remaining.set_selected(s,-1.0)
  return map_data_remaining

def get_lower(lower_bounds,lower):
  new_lower=[]
  for i in xrange(3):
    if lower_bounds[i] is None:
      new_lower.append(lower[i])
    elif lower[i] is None:
      new_lower.append(lower_bounds[i])
    else:
      new_lower.append(min(lower_bounds[i],lower[i]))
  return new_lower

def get_upper(upper_bounds,upper):
  new_upper=[]
  for i in xrange(3):
    if upper_bounds[i] is None:
      new_upper.append(upper[i])
    elif upper[i] is None:
      new_upper.append(upper_bounds[i])
    else:
      new_upper.append(max(upper_bounds[i],upper[i]))
  return new_upper

def get_bounds(min_b=None,max_b=None,ncs_group_obj=None,id=None):
  orig_id=ncs_group_obj.original_id_from_id[id]
  lower=min_b[orig_id]
  upper=max_b[orig_id]
  return lower,upper

def get_selected_and_related_regions(params,
    ncs_group_obj=None):
  # Identify all points in the targeted regions
  bool_selected_regions=get_bool_mask_of_regions(ncs_group_obj=ncs_group_obj,
     region_list=ncs_group_obj.selected_regions,
     expand_size=params.segmentation.expand_size)
  # and all points in NCS-related copies (to be excluded)
  bool_ncs_related_mask=get_bool_mask_of_regions(ncs_group_obj=ncs_group_obj,
       region_list=ncs_group_obj.ncs_related_regions)

  min_b,max_b=ncs_group_obj.co.get_blobs_boundaries_tuples()
  lower_bounds=[None,None,None]
  upper_bounds=[None,None,None]
  for id in ncs_group_obj.selected_regions:
    lower,upper=get_bounds(min_b=min_b,max_b=max_b,ncs_group_obj=ncs_group_obj,
       id=id)
    lower_bounds=get_lower(lower_bounds,lower)
    upper_bounds=get_upper(upper_bounds,upper)

  return bool_selected_regions,bool_ncs_related_mask,lower_bounds,upper_bounds

def adjust_bounds(params,lower_bounds,upper_bounds,map_data=None,out=sys.stdout):
  # range is lower_bounds to upper_bounds
  lower_bounds=list(lower_bounds)
  upper_bounds=list(upper_bounds)
  if params.output_files.box_buffer is None: params.output_files.box_buffer=0
  for i in xrange(3):
    if lower_bounds[i] is None: lower_bounds[i]=0
    if upper_bounds[i] is None: upper_bounds[i]=0
    lower_bounds[i]-=params.output_files.box_buffer
    lower_bounds[i]=max(0,lower_bounds[i])
    upper_bounds[i]+=params.output_files.box_buffer
    upper_bounds[i]=min(map_data.all()[i],upper_bounds[i])


  """
  print >>out,"\nRange:  X:(%6d,%6d)    Y:(%6d,%6d)    Z:(%6d,%6d)" %(
     lower_bounds[0],upper_bounds[0],
     lower_bounds[1],upper_bounds[1],
     lower_bounds[2],upper_bounds[2])
  """

  return lower_bounds,upper_bounds

def write_region_maps(params,
    ncs_group_obj=None,
    map_data=None,
    crystal_symmetry=None,
    remainder_ncs_group_obj=None,
    regions_to_skip=None,
    out=sys.stdout):
  remainder_regions_written=[]
  map_files_written=[]
  if not ncs_group_obj:
    return map_files_written,remainder_regions_written

  min_b,max_b=ncs_group_obj.co.get_blobs_boundaries_tuples()
  if remainder_ncs_group_obj:
    r_min_b,r_max_b=remainder_ncs_group_obj.co.get_blobs_boundaries_tuples()


  for id in ncs_group_obj.selected_regions:
    if regions_to_skip and id in regions_to_skip:
      print >>out,"Skipping remainder region %d (already written out)" %(id)
      continue
    print >>out,"Writing region %d" %(id),
    bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand=ncs_group_obj.original_id_from_id[id],
        expand_size=params.segmentation.expand_size)

    s = (bool_region_mask==True)

    lower_bounds,upper_bounds=get_bounds(
      min_b=min_b,max_b=max_b,ncs_group_obj=ncs_group_obj,id=id)

    if remainder_ncs_group_obj:
      for remainder_id in remainder_ncs_group_obj.remainder_id_dict.keys():
        if remainder_ncs_group_obj.remainder_id_dict[remainder_id]==id:
          remainder_regions_written.append(remainder_id)
          print >>out,"(including remainder region %d)" %(remainder_id),
          remainder_bool_region_mask = remainder_ncs_group_obj.co.expand_mask(
           id_to_expand=remainder_ncs_group_obj.original_id_from_id[remainder_id],
           expand_size=params.segmentation.expand_size)
          s|= (remainder_bool_region_mask==True)
          lower,upper=get_bounds(
            min_b=r_min_b,max_b=r_max_b,ncs_group_obj=remainder_ncs_group_obj,
            id=remainder_id)
          lower_bounds=get_lower(lower_bounds,lower)
          upper_bounds=get_upper(upper_bounds,upper)

    region_mask = map_data.deep_copy()
    region_mask = region_mask.set_selected(s,1)
    region_mask = region_mask.set_selected(~s,0)
    local_map_data=map_data.deep_copy()
    local_map_data=local_map_data * region_mask.as_double()

    # Now cut down the map to the size we want
    lower_bounds,upper_bounds=adjust_bounds(params,lower_bounds,upper_bounds,
      map_data=map_data,out=out)
    box_map,box_crystal_symmetry=cut_out_map(
       map_data=local_map_data, crystal_symmetry=crystal_symmetry,
       min_point=lower_bounds, max_point=upper_bounds)

    if remainder_ncs_group_obj:
      text=""
    else:
      text="_r"
    base_file='map%s_%d.ccp4' %(text, id)
    if params.output_files.output_map_directory:
      file_name=os.path.join(params.output_files.output_map_directory,base_file)
    else:
      file_name=base_file
    write_ccp4_map(box_crystal_symmetry,file_name, box_map)
    print >>out,"to %s" %(file_name)
    map_files_written.append(file_name)
  return map_files_written,remainder_regions_written

def write_output_files(params,
    crystal_symmetry=None,
    map_data=None,
    map_data_remaining=None,
    ncs_group_obj=None,
    remainder_ncs_group_obj=None,
    pdb_hierarchy=None,
    out=sys.stdout):

  # Write out mask and map representing one NCS copy and none of
  #   other NCS copies.  Expand the mask to include neighboring points (but
  #   not those explicitly in other NCS copies

  bool_selected_regions,bool_ncs_related_mask,lower_bounds,upper_bounds=\
     get_selected_and_related_regions(
      params,ncs_group_obj=ncs_group_obj)
  s=  (bool_ncs_related_mask==True)

  # Add in remainder regions if present
  if remainder_ncs_group_obj:
    bool_remainder_selected_regions,bool_remainder_ncs_related_mask,\
      remainder_lower_bounds,remainder_upper_bounds=\
       get_selected_and_related_regions(
       params,ncs_group_obj=remainder_ncs_group_obj)

    lower_bounds=get_lower(lower_bounds,remainder_lower_bounds)
    upper_bounds=get_upper(upper_bounds,remainder_upper_bounds)

    s_remainder_au =  (bool_remainder_selected_regions==True)
    bool_selected_regions=bool_selected_regions.set_selected(
       s_remainder_au,True)

    s |=  (bool_remainder_ncs_related_mask==True)

  # Now create NCS mask by eliminating all points in target (expanded) in
  #   NCS-related copies

  bool_selected_regions=bool_selected_regions.set_selected(s,False)

  lower_bounds,upper_bounds=adjust_bounds(params,lower_bounds,upper_bounds,
    map_data=map_data,out=out)

  print >>out,\
     "\nMaking two types of maps for AU of NCS mask and map with "+\
      "buffer of %d grid units \nin each direction around AU" %(
      params.output_files.box_buffer)
  print >>out,"Both types of maps have the same origin and overlay on %s" %(
   params.output_files.shifted_map_file) 
      
  print >>out,\
     "\nThe standard maps (%s, %s) have the \noriginal cell dimensions." %(
   params.output_files.au_mask_output_file,
   params.output_files.au_map_output_file)+\
   "\nThese maps show only the unique (NCS AU) part of the map."

  print >>out,\
     "\nThe cut out box_maps (%s, %s) have \nsmaller cell dimensions." %(
      params.output_files.box_mask_file,
      params.output_files.box_map_file,) +\
   "\nThese maps also show only the unique part of the map and have this"+\
   "\nunique part cut out.\n"
 

  # Write out mask and map of NCS AU with shifted origin but 
  #  initial crystal_symmetry

  mask_data_ncs_au=get_bool_mask_as_int(
     ncs_group_obj=ncs_group_obj,mask_as_bool=bool_selected_regions)

  if params.output_files.au_mask_output_file: # Write out the mask (as int)
    write_ccp4_map(crystal_symmetry,
      params.output_files.au_mask_output_file,mask_data_ncs_au)
    print >>out,\
       "Output NCS AU mask:  %s" %(
      params.output_files.au_mask_output_file)

  # Now write out NCS AU mask and map cut out, with box_crystal_symmetry

  map_data_ncs_au=map_data.deep_copy()
  s=(bool_selected_regions==True)
  map_data_ncs_au=map_data_ncs_au.set_selected(~s,-1.0)

  if params.output_files.au_map_output_file: # Write out the NCS au of density
    write_ccp4_map(crystal_symmetry,params.output_files.au_map_output_file,
      map_data_ncs_au)
    print >>out,\
       "Output NCS AU map:  %s" %(
      params.output_files.au_map_output_file)

  box_mask_ncs_au,box_crystal_symmetry=cut_out_map(
       map_data=mask_data_ncs_au.as_double(), crystal_symmetry=crystal_symmetry,
       min_point=lower_bounds, max_point=upper_bounds)

  if params.output_files.box_mask_file:
    # write out NCS mask and mask as box_mask (cut out the region they enclose)
    write_ccp4_map(
     box_crystal_symmetry,params.output_files.box_mask_file,
      box_mask_ncs_au)
    print >>out,\
      "Output NCS au as box (cut out) mask:  %s " %(
      params.output_files.box_mask_file)

  if params.output_files.box_map_file:
    # write out NCS map and mask as box_map (cut out the region they enclose)
    box_map_ncs_au,box_crystal_symmetry=cut_out_map(
       map_data=map_data_ncs_au.as_double(), crystal_symmetry=crystal_symmetry,
       min_point=lower_bounds, max_point=upper_bounds)
    write_ccp4_map(box_crystal_symmetry,params.output_files.box_map_file,
      box_map_ncs_au)
    print >>out,\
      "Output NCS au as box (cut out) map:  %s " %(
      params.output_files.box_map_file)

  # Write out all the selected regions
  print >>out,"\nWriting out region maps. "+\
    "These superimpose on the NCS AU map \nand "+\
    "mask %s,%s\n" %(
      params.output_files.box_map_file,
      params.output_files.box_mask_file,)

  map_files_written,remainder_regions_written=write_region_maps(params,
    map_data=map_data,
    crystal_symmetry=crystal_symmetry,
    ncs_group_obj=ncs_group_obj,
    remainder_ncs_group_obj=remainder_ncs_group_obj,
    out=out)

  # and pick up the remainder regions not already written
  remainder_map_files_written,dummy_remainder=write_region_maps(params,
    map_data=map_data,
    crystal_symmetry=crystal_symmetry,
    ncs_group_obj=remainder_ncs_group_obj,
    regions_to_skip=remainder_regions_written,
    out=out)
  map_files_written+=remainder_map_files_written
  return map_files_written

def write_intermediate_maps(params,
    crystal_symmetry=None,
    map_data=None,
    map_data_remaining=None,
    ncs_group_obj=None,
    out=sys.stdout):

  if map_data_remaining and params.output_files.remainder_map_file:
    write_ccp4_map(crystal_symmetry,params.output_files.remainder_map_file,
      map_data_remaining)
    print >>out,"Wrote output remainder map to %s" %(
       params.output_files.remainder_map_file)

  if params.segmentation.write_all_regions:
    for id in ncs_group_obj.selected_regions:
      region_mask=ncs_group_obj.edited_mask.deep_copy()
      s = (ncs_group_obj.edited_mask == -1)
      s |= (ncs_group_obj.edited_mask==id)
      region_mask = region_mask.set_selected(s,1)
      region_mask = region_mask.set_selected(~s,0)

      write_ccp4_map(crystal_symmetry,'mask_%d.ccp4' %id, region_mask)
      print >>out,"Wrote output mask for region %d to %s" %(id,
        "mask_%d.ccp4" %(id))


def iterate_search(params,
      map_data_remaining=None,
      map_data=None,
      n_residues=None,
      solvent_fraction=None,
      ncs_obj=None,
      ncs_group_obj=None,
      scattered_points=None,
      crystal_symmetry=None,
      out=sys.stdout):


  # Write out intermediate maps if desired
  if params.output_files.write_intermediate_maps:
    write_intermediate_maps(params,
      crystal_symmetry=crystal_symmetry,
      map_data=map_data,
      map_data_remaining=map_data_remaining,
      ncs_group_obj=ncs_group_obj,
      out=out)

  from copy import deepcopy
  new_params=deepcopy(params)
  new_params.segmentation.iterate_with_remainder=False
  new_params.segmentation.density_threshold=None
  new_params.output_files.write_output_maps=False
  new_params.output_files.output_info_file=None
  if params.output_files.write_intermediate_maps:
    new_params.output_files.au_map_output_file=\
     params.output_files.au_map_output_file[:-5]+"_cycle_2.ccp4"
    new_params.output_files.au_mask_output_file=\
      params.output_files.au_mask_output_file[:-5]+"_cycle_2.ccp4"
  else:
    new_params.output_files.au_map_output_file=None
    new_params.output_files.au_mask_output_file=None

  fraction=0.2
  new_n_residues=int(n_residues*fraction)
  new_solvent_fraction=1- (1-solvent_fraction)*fraction

  print >>out,"\nIterating with remainder density"
  # NOTE: do not include pdb_hierarchy here unless you deep_copy it
  remainder_ncs_group_obj,dummy_remainder=run(None,params=new_params,
    map_data=map_data_remaining,
    ncs_obj=ncs_obj,
    target_scattered_points=scattered_points,
    solvent_fraction=new_solvent_fraction,
    n_residues=new_n_residues,
    crystal_symmetry=crystal_symmetry,
    out=out)
  if not remainder_ncs_group_obj: # Nothing to do
    return None

  # Combine the results to get remainder_id_dict
  #   remainder_id_dict[id_remainder]=id_nearby

  remainder_ncs_group_obj=combine_with_iteration(params,
     map_data=map_data,
     crystal_symmetry=crystal_symmetry,
     ncs_group_obj=ncs_group_obj,
     remainder_ncs_group_obj=remainder_ncs_group_obj,
     out=out)

  return remainder_ncs_group_obj

def combine_with_iteration(params,
    map_data=None,
    crystal_symmetry=None,
    ncs_group_obj=None,
    remainder_ncs_group_obj=None,
    out=sys.stdout):

  if not ncs_group_obj.selected_regions or \
      not remainder_ncs_group_obj or not remainder_ncs_group_obj.selected_regions:
    return None

  # see if any regions in ncs_obj overlap with remainder_ncs_group_obj...
  #   if so, combine
  remainder_id_dict={}
  for id_remainder in remainder_ncs_group_obj.selected_regions:
    best_id=None
    best_overlaps=None
    for id in ncs_group_obj.selected_regions:
      bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand=ncs_group_obj.original_id_from_id[id],
        expand_size=params.segmentation.expand_size+1) # just touching
      s = (bool_region_mask== True)
      s &=  (remainder_ncs_group_obj.edited_mask==id_remainder)
      overlaps=s.count(True)
      if best_overlaps is None or overlaps>best_overlaps:
        best_overlaps=overlaps
        best_id=id
    if best_overlaps:
      print >>out,\
        "\nCombining remainder id %d with original id %d (overlaps=%d)" %(
        id_remainder,best_id,best_overlaps)
      remainder_id_dict[id_remainder]=best_id
  remainder_ncs_group_obj.remainder_id_dict=remainder_id_dict
  return remainder_ncs_group_obj


def cut_out_map(map_data=None, crystal_symmetry=None,
    min_point=None,max_point=None):
  from cctbx import uctbx
  from cctbx import maptbx
  na = map_data.all() # tuple with dimensions!!!
  for i in range(3):
    assert min_point[i] >= 0
    assert max_point[i] <= na[i]
  new_map_data = maptbx.copy(map_data, tuple(min_point), tuple(max_point))
  # shrink unit cell, angles are the same
  shrunk_uc = []
  for i in range(3):
    shrunk_uc.append(
     crystal_symmetry.unit_cell().parameters()[i] * new_map_data.all()[i]/na[i] )
  uc_params=crystal_symmetry.unit_cell().parameters()
  new_unit_cell_box = uctbx.unit_cell(
    parameters=(shrunk_uc[0],shrunk_uc[1],shrunk_uc[2],
        uc_params[3],uc_params[4],uc_params[5]))
  new_crystal_symmetry=crystal.symmetry(
    unit_cell=new_unit_cell_box,space_group=crystal_symmetry.space_group())
  return new_map_data, new_crystal_symmetry

def apply_shift_to_pdb_hierarchy(
    origin_shift=None,
    crystal_symmetry=None,
    pdb_hierarchy=None,out=sys.stdout):

  if origin_shift is not None:
    sites_cart=pdb_hierarchy.atoms().extract_xyz()
    sites_cart_shifted=sites_cart+\
      flex.vec3_double(sites_cart.size(), origin_shift)
    pdb_hierarchy.atoms().set_xyz(sites_cart_shifted)

  return pdb_hierarchy

def apply_origin_shift(origin_shift=None,
    ncs_object=None,
    pdb_hierarchy=None,
    crystal_symmetry=None,
    map_data=None,
    shifted_map_file=None,
    shifted_pdb_file=None,
    shifted_ncs_file=None,
    min_point=None, # ZZZ
    max_point=None,
    out=sys.stdout):

  if origin_shift: # Note origin shift does not change crystal_symmetry
    if shifted_map_file:
      write_ccp4_map(crystal_symmetry,shifted_map_file,
      map_data)
      print >>out,"Wrote shifted map to %s" %(
        shifted_map_file)

    if pdb_hierarchy:
      pdb_hierarchy=apply_shift_to_pdb_hierarchy(
       origin_shift=origin_shift,
       crystal_symmetry=crystal_symmetry,
       pdb_hierarchy=pdb_hierarchy,
       out=out)
      
    if shifted_pdb_file: 
      import iotbx.pdb
      f=open(shifted_pdb_file,'w')
      print >>f, iotbx.pdb.format_cryst1_record(
         crystal_symmetry=crystal_symmetry)
      print >>f,pdb_hierarchy.as_pdb_string()
      f.close() 
      print >>out,"Wrote shifted pdb file to %s" %(
        shifted_pdb_file)

    from scitbx.math import  matrix
    ncs_object=ncs_object.coordinate_offset(
       coordinate_offset=matrix.col(origin_shift))
    #ncs_object.display_all(log=out)
    if shifted_ncs_file:
      ncs_object.format_all_for_group_specification(
         file_name=shifted_ncs_file)
      print >>out,"Wrote NCS operators for shifted map to %s" %(
       shifted_ncs_file)
    
  return ncs_object,pdb_hierarchy

def run(args,
     params=None,
     map_data=None,
     ncs_obj=None,
     solvent_fraction=None,
     n_residues=None,
     crystal_symmetry=None,
     target_scattered_points=None,
     out=sys.stdout):

  if params and  map_data and ncs_obj and solvent_fraction and \
      n_residues and crystal_symmetry:  # have things ready to go
    origin_shift=None
  else:
    # get the parameters
    params,crystal_symmetry,map_data,origin_shift,pdb_hierarchy=get_params(
        args,out=out)

    # read and write the ncs (Normally point-group NCS)
    ncs_obj=get_ncs(params)

    if origin_shift:
      print >>out,"\nShifting map, model and NCS based on origin shift"
      ncs_obj,pdb_hierarchy=apply_origin_shift(
        shifted_map_file=params.output_files.shifted_map_file,
        shifted_ncs_file=params.output_files.shifted_ncs_file,
        shifted_pdb_file=params.output_files.shifted_pdb_file,
        origin_shift=origin_shift,
        ncs_object=ncs_obj,
        pdb_hierarchy=pdb_hierarchy,
        crystal_symmetry=crystal_symmetry,
        map_data=map_data,
        out=out)
    else:
      shifted_map_file=params.input_files.ccp4_map_file # but do not overwrite

    # get the chain types and therefore (using ncs_copies) volume fraction
    solvent_fraction,n_residues=get_solvent_fraction(
       params,crystal_symmetry=crystal_symmetry,
       ncs_object=ncs_obj,out=out)


  # get connectivity  (conn=connectivity_object.result)
  co,conn,sorted_by_volume,min_b,max_b=get_connectivity(params,
     map_data=map_data,
     ncs_object=ncs_obj,
     n_residues=n_residues,
     ncs_copies=ncs_obj.max_operators(),
     solvent_fraction=solvent_fraction,out=out)

  # Check to see which regions are in more than one au of the NCS
  #   and set them aside.  Group ncs-related regions together

  ncs_group_obj=identify_ncs_regions(
     params,sorted_by_volume=sorted_by_volume,
     unit_cell=crystal_symmetry.unit_cell(),
     crystal_symmetry=crystal_symmetry,
     co=co,
     conn_obj=conn,
     ncs_obj=ncs_obj,
     origin_shift=origin_shift,
     out=out)
  if not ncs_group_obj:  # nothing to do
    return None,None

  # Choose one region or group of regions from each ncs_group in the list
  #  Optimize the closeness of centers

  # Select group of regions that are close together and represent one
  ncs_group_obj,scattered_points=\
     select_regions_in_au(
     params,
     ncs_group_obj=ncs_group_obj,
     target_scattered_points=target_scattered_points,
     out=out)

  # write out mask and map for all the selected regions...

  map_data_remaining=create_remaining_mask_and_map(params,
    ncs_group_obj=ncs_group_obj,
    map_data=map_data,
    crystal_symmetry=crystal_symmetry,
    out=out)

  # Iterate if desired
  if params.segmentation.iterate_with_remainder:
    remainder_ncs_group_obj=iterate_search(params,
      map_data=map_data,
      map_data_remaining=map_data_remaining,
      n_residues=n_residues,
      solvent_fraction=solvent_fraction,
      ncs_obj=ncs_obj,
      ncs_group_obj=ncs_group_obj,
      scattered_points=scattered_points,
      crystal_symmetry=crystal_symmetry,
      out=out)
  else:
    remainder_ncs_group_obj=None

  # Write out final maps
  if params.output_files.write_output_maps:
    map_files_written=write_output_files(params,
      crystal_symmetry=crystal_symmetry,
      map_data=map_data,
      map_data_remaining=map_data_remaining,
      ncs_group_obj=ncs_group_obj,
      remainder_ncs_group_obj=remainder_ncs_group_obj,
      pdb_hierarchy=pdb_hierarchy,
      out=out)
    ncs_group_obj.set_map_files_written(map_files_written)
  else:
    map_files_written=[]

  if params.output_files.output_info_file:
    from libtbx import easy_pickle
    info_obj=ncs_group_obj.as_info_object()
    easy_pickle.dump( params.output_files.output_info_file,
       info_obj)

  return ncs_group_obj,remainder_ncs_group_obj


if __name__=="__main__":
  run(args=sys.argv[1:])
