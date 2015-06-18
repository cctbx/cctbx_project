from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.place_model_in_map
import sys,os
import iotbx.phil
from libtbx.utils import Sorry,null_out

# Take a PDB file and rotate/translate it to fit into a .ccp4 map
# Quick by using low-resolution MR. Can make it more accurate with
# high-resolution MR or rigid-body refinement

# First check to see if simple translation will do it...

master_phil = iotbx.phil.parse("""
  input_files {

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    pdb_file = None
      .type = path
      .help = Input PDB file to match to the map
      .short_caption = Input PDB file

  }
  directories {
    temp_dir = temp_place_model_in_map
      .type = path
      .help = temporary directory
      .short_caption = temporary directory
  }
  output_files {
    pdb_out = placed_model.pdb
      .type = path
      .help = Output PDB file matched to map
      .short_caption = Output PDB file
  }
  place_model {
     resolution = None
       .type = float
       .help = Resolution of the map.  Required input.
       .short_caption = Map resolution
       .style = resolution
     molecular_mass = None
       .type = float
       .help = Molecular mass of molecule.  Normally estimated from your model.
       .short_caption = Molecular mass

     solvent_content = 0.85
       .type = float
       .help = Solvent fraction of the cell.  If this is density cut out from \
          a bigger cell, you can specify the fraction of the volume of \
          this cell that is taken up by the macromolecule. \
          Values go from 0 to 1.
       .short_caption = Solvent content
     resolution_factor = 2
       .type = float
       .help = Resolution for finding placement of model will be resolution \
          times the resolution_factor
       .short_caption = Resolution factor
  }
  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output
      clean_up = True
        .type = bool
        .help = Remove temp_dir when finished
        .short_caption = Remove temp_dir when finished
  }
""", process_includes=True)

master_params = master_phil

def get_params(args,out=sys.stdout):
  command_line = iotbx.phil.process_command_line_with_files(
    pdb_file_def="input_files.pdb_file",
    map_file_def="input_files.map_file",
    args=args,
    master_phil=master_phil)
  params = command_line.work.extract()

  # make sure we have resolution
  if params.place_model.resolution is None:
    for arg in args:
      if arg.find("=")<0 and not os.path.isfile(arg):
        try:
          params.place_model.resolution=float(arg)
        except Exception: pass  # no need

  if not params.place_model.resolution:
    raise Sorry("Need to specify resolution=xxx for place_model_in_map")

  print >>out,"\nPlace_model_in_map...Rotate/translate a model to match a map"
  master_phil.format(python_object=params).show(out=out)

  return params

def guess_mw(file_name):
  # Guess MW from the PDB file
  import iotbx.pdb
  from cctbx.array_family import flex
  text=open(file_name).read()
  if not text:
    raise Sorry(
      "Nothing found in the PDB file %s" %(params.input_files.pdb_file))
  pdb_inp=iotbx.pdb.input(source_info="",
       lines=flex.split_lines(text))
  hierarchy=pdb_inp.construct_hierarchy()
  asc=hierarchy.atom_selection_cache()
  atom_selection="not element H"
  sel = asc.selection(string = atom_selection)
  hierarchy = hierarchy.select(sel)
  n_atoms=hierarchy.overall_counts().n_atoms
  mw=n_atoms*6.7 # Hendrickson, W.A. (2014) Quarterly Rev. Biophys. 47, 49-93.
  print "\nTotal of %d atoms. Molecular mass is about %.0f Da"%(
     hierarchy.overall_counts().n_atoms,mw)
  return hierarchy,mw,n_atoms

def run(args,out=sys.stdout):

  # get the parameters
  params=get_params(args,out=out)

  # make a temporary directory
  if not os.path.isdir(params.directories.temp_dir):
    os.mkdir(params.directories.temp_dir)

  # choose resolution for FFT
  d_min=params.place_model.resolution*params.place_model.resolution_factor
  print >>out,"\nWorking resolution will be %7.1f A" %(d_min)

  # Get the map as structure factors and write to a temp file
  mtz_file_name=os.path.join(params.directories.temp_dir,
    "map_as_structure_factors.mtz")
  from mmtbx.command_line.map_to_structure_factors import run as mtsf
  mtsf_args=[params.input_files.map_file,
    "output_file_name=%s" %(mtz_file_name),
    "d_min=%s" %(str(d_min))]
  mtsf(mtsf_args,nohl=True,out=out)

  hierarchy,mw,n_atoms=guess_mw(params.input_files.pdb_file)
  if params.place_model.molecular_mass is None:
    params.place_model.molecular_mass=mw

  # write out the file again so Phaser can read it (phaser can't read cif)
  f=open(os.path.join(params.directories.temp_dir,'starting_model.pdb'),'w')
  pdb_for_phaser=f.name
  print >>f,hierarchy.as_pdb_string()
  f.close()
  print >>out,"\nStarting model copied to %s" %(pdb_for_phaser)
  # Now just do MR using the data in our temporary mtz file
  root="mr"
  phaser_file_name=root+".1.pdb"
  phaser_args=[
   "mode=MR_AUTO",
   "%s" %(mtz_file_name),
   "labin=F-obs,SIGF-obs",
   "%s" %(pdb_for_phaser),
   "mol_weight=%.0f" %(params.place_model.molecular_mass),
   "model_identity=90",
   "component_copies=1",
   "search.copies=1",
   #"output_dir=%s" %(params.directories.temp_dir),  # does not work
   "keywords.general.root=%s" %(root),
   "keywords.general.xyzout=True",
   "keywords.general.hklout=False",
   "sgalternative.select=None",
   ]
  from phaser.phenix_interface import driver
  if params.control.verbose:
    print >>out,"Phaser args: %s" %(" ".join(phaser_args))
  if params.control.verbose:
    local_out=out
  else:
    local_out=null_out()
  print >>out,"\nRunning MR to find orientation of model..."
  driver.run(phaser_args,out=local_out)
  if os.path.isfile(phaser_file_name):
    print >>out,"MR solution (not offset yet) is in %s" %(phaser_file_name)
  else:
    raise Sorry("Unable to place model")

  # Now move the model (translate) to match map (the MR was in P1)
  from phenix.command_line.get_cc_mtz_pdb import get_cc_mtz_pdb
  cc_args=[
   "pdb_in=%s" %(phaser_file_name),
   "mtz_in=%s" %(mtz_file_name),
   "offset_pdb=%s" %(params.output_files.pdb_out),
   "quick=True",
   "resolution=%7.1f" %(d_min),
   "temp_dir=%s" %(params.directories.temp_dir),
   "output_dir=%s" %(os.getcwd()),
   ]
  print >>out,"\nTranslating MR model to match map..."
  cc_mtz_pdb=get_cc_mtz_pdb(cc_args,copy_extra_columns=False,out=local_out)
  cc_value=cc_mtz_pdb.found_region

  assert os.path.isfile(params.output_files.pdb_out)
  print >>out,"\nModel matching map is in %s\n" %(params.output_files.pdb_out)
  print >>out,"CC of model to map is %5.2f" %(cc_value)

  if params.control.clean_up:
    print >>out,"Removing temporary directory %s" %(params.directories.temp_dir)
    from libtbx.clear_paths import \
       remove_or_rename_files_and_directories_if_possible
    remove_or_rename_files_and_directories_if_possible(
       paths=[params.directories.temp_dir])

  return params.output_files.pdb_out,cc_value
if __name__=="__main__":
  args=sys.argv[1:]
  run(args)
