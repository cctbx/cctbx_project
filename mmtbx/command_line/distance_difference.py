# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# TODO need a good test for this

import libtbx.phil.command_line
from libtbx.utils import Sorry, Usage, null_out
import os
import sys

calculate_matrix_phil = """
model_1 = None
  .type = path
  .optional = False
chain_1 = None
  .type = str
model_2 = None
  .type = path
  .optional = False
chain_2 = None
  .type = str
"""

master_phil = libtbx.phil.parse("""
calculate_matrix
  .caption = This tool calculates a distance-difference matrix for two \
    protein chains.  The resulting plot shows the change in interatomic \
    distances (between C-alpha atoms), illustrating correlated motions and \
    rigid bodies in the structures.  The input structures are assumed to be \
    identical in sequence numbering, but point mutations and missing regions \
    are allowed.
  .short_caption = Calculate distance-difference matrix
  .style = auto_align box caption_img:icons/custom/distance_difference.png
{
%s
display_plot = False
  .type = bool
  .style = hidden
}""" % calculate_matrix_phil)

def usage () :
  raise Usage("""
mmtbx.distance_difference model1.pdb CHAIN_1 model2.pdb CHAIN_2
mmtbx.distance_difference model_1=... model_2=... chain_1=... chain_2=...

other options:
  --display_plot or display_plot=True : display plot interactively

Calculate and plot a distance difference matrix for two (near-identical)
protein chains.  Chain IDs are optional for models that have only one protein
chain.
""")

def run (args=(), params=None, out=None, display_plot=False) :
  if (out is None) :
    out = sys.stdout
  user_phil = []
  pdb_files = []
  chain_ids = []
  interp = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="calculate_matrix")
  if (params is None) :
    if (len(args) == 0) :
      usage()
    import iotbx.pdb
    for arg in args :
      if os.path.isfile(arg) :
        if (iotbx.pdb.is_pdb_file(arg)) :
          pdb_files.append(arg)
        else :
          user_phil.append(libtbx.phil.parse(file_name=arg))
      elif (len(arg) <= 2) :
        chain_ids.append(arg)
      else :
        if (arg.startswith("--")) :
          arg = arg[2:] + "=True"
        try :
          arg_phil = interp.process_arg(arg)
        except RuntimeError :
          pass
        else :
          user_phil.append(arg_phil)
          continue
    for i, file_name in enumerate(pdb_files) :
      i += 1
      user_phil.append(interp.process_arg("model_%d=\"%s\"" % (i, file_name)))
    for i, chain_id in enumerate(chain_ids) :
      i += 1
      user_phil.append(interp.process_arg("chain_%d=\"%s\"" % (i, chain_id)))
    working_phil = master_phil.fetch(sources=user_phil)
    params = working_phil.extract()
    validate_params(params)
  params = params.calculate_matrix
  ddm = calculate_matrix(params, log=out)
  title = "Distance-difference matrix: %s:%s vs. %s:%s" % (
    os.path.basename(params.model_1), params.chain_1,
    os.path.basename(params.model_2), params.chain_2)
  ddm.set_title(title)
  if (display_plot) and (params.display_plot) :
    try :
      display_plot_pylab(ddm)
    except ImportError :
      raise Sorry("matplotlib is not installed - can't generate plot image.")
    except Exception, e :
      print >> out, "Oops!  Can't display an interactive plot:"
      print >> out, "  %s" % str(e)
  elif (display_plot) : # only if not run from GUI!
    try :
      import matplotlib
    except ImportError :
      raise Sorry("matplotlib is not installed - can't generate plot image.")
    matplotlib.use('Agg')
    display_plot_pylab(ddm, savefig=True)
    print "Saved plot as distance_difference.png"
  return ddm

def calculate_matrix (params, log=None) :
  if (log is None) : log = null_out()
  from iotbx import file_reader
  pdb_1 = file_reader.any_file(params.model_1, force_type="pdb")
  pdb_1.check_file_type("pdb")
  hierarchy_1 = pdb_1.file_object.construct_hierarchy()
  hierarchy_1.atoms().reset_i_seq()
  pdb_2 = file_reader.any_file(params.model_2, force_type="pdb")
  pdb_2.check_file_type("pdb")
  hierarchy_2 = pdb_2.file_object.construct_hierarchy()
  hierarchy_2.atoms().reset_i_seq()
  for k, hierarchy in enumerate([hierarchy_1,hierarchy_2]) :
    k += 1
    if (getattr(params, "chain_%d" % k) is None) :
      try :
        chain_id = find_single_protein_chain(hierarchy)
        if (chain_id is None) :
          raise Sorry(("The file %s does not appear to have any protein "+
            "chains!  If you think this is an error, you should explicitly "+
            "specify the parameter chain_%d.  (Nucleic acids are not yet "+
            "supported, sorry.)") % (getattr(params, "model_%d" % k), k))
        setattr(params, "chain_%d" % k, chain_id)
      except RuntimeError :
        raise Sorry(("The file %s has more than one protein chain, and no "+
          "explicit chain ID to use was specified.") %
            getattr(params, "model_%d" % k))
      except AssertionError :
        raise Sorry(("The file %s contains multiple models.  You can use "+
          "iotbx.pdb.split_models to extract the individual models from the "+
          "file (also available in the PHENIX GUI).") %
            getattr(params, "model_%d" % k))
  return distance_difference_matrix(
    hierarchy_1=hierarchy_1,
    chain_id_1=params.chain_1,
    hierarchy_2=hierarchy_2,
    chain_id_2=params.chain_2,
    log=log)

# TODO any polymer, not just protein
def find_single_protein_chain (hierarchy) :
  assert (len(hierarchy.models()) == 1)
  chain_id = None
  for chain in hierarchy.models()[0].chains() :
    main_conf = chain.conformers()[0]
    if (main_conf.is_protein()) :
      if (chain_id is not None) :
        raise RuntimeError("More than one protein chain in hierarchy.")
      chain_id = chain.id
  return chain_id

class distance_difference_matrix (object) :
  def __init__ (self,
                hierarchy_1,
                chain_id_1,
                hierarchy_2,
                chain_id_2,
                log=None) :
    assert (not None in [hierarchy_1, chain_id_1, hierarchy_2, chain_id_2])
    assert (len(hierarchy_1.models())==1) and (len(hierarchy_2.models())==1)
    assert (not hierarchy_1.atoms().extract_i_seq().all_eq(0))
    assert (not hierarchy_2.atoms().extract_i_seq().all_eq(0))
    from scitbx.array_family import flex
    if (log is None) :
      log = null_out()
    self.matching_resids = []
    chain_1 = chain_2 = None
    for chain in hierarchy_1.models()[0].chains() :
      if (chain.id == chain_id_1) :
        main_conf = chain.conformers()[0]
        if (main_conf.is_protein()) :
          if (chain_1 is not None) :
            raise Sorry("Multiple protein chains with ID '%s' in hierarchy_1")
          else :
            chain_1 = chain
    for chain in hierarchy_2.models()[0].chains() :
      if (chain.id == chain_id_2) :
        main_conf = chain.conformers()[0]
        if (main_conf.is_protein()) :
          if (chain_2 is not None) :
            raise Sorry("Multiple protein chains with ID '%s' in hierarchy_1")
          else :
            chain_2 = chain
    if (chain_1 is None) :
      raise Sorry("Can't find protein chain '%s' in 1st model." % chain_id_1)
    if (chain_2 is None) :
      raise Sorry("Can't find protein chain '%s' in 2nd model." % chain_id_2)
    resids = [ [], [] ]
    selections = [ flex.size_t(), flex.size_t() ]
    for i, chain in enumerate([chain_1, chain_2]) :
      for residue_group in chain.residue_groups() :
        atom_groups = residue_group.atom_groups()
        resid = residue_group.resid()
        if (len(atom_groups) > 1) :
          print >> log, "Warning: multiple conformers in model %d at %s" % (i+1,
            resid)
        for atom in atom_groups[0].atoms() :
          if (atom.name == " CA ") :
            resids[i].append(resid)
            selections[i].append(atom.i_seq)
            break
    for i,j in [ (0,1), (1,0) ] :
      k = 0
      assert (len(resids[i]) == len(selections[i]))
      while (k < len(resids[i])) :
        if (not resids[i][k] in resids[j]) :
          del resids[i][k]
          del selections[i][k]
        else :
          k += 1
    assert (len(resids[0]) == len(resids[1]))
    assert (len(selections[0]) == len(selections[1]))
    if (len(resids[0]) == 0) :
      raise Sorry("No matching residues in model 1!  Make sure it is "+
        "really a protein chain with CA atoms for each residue.")
    if (len(resids[1]) == 0) :
      raise Sorry("No matching residues in model 2!  Make sure it is "+
        "really a protein chain with CA atoms for each residue.")
    sites_1 = hierarchy_1.atoms().extract_xyz()
    sites_1 = sites_1.select(selections[0])
    sites_2 = hierarchy_2.atoms().extract_xyz()
    sites_2 = sites_2.select(selections[1])
    import scitbx.math
    self.m = scitbx.math.distance_difference_matrix(sites_1, sites_2)
    self.resids_1 = resids[0]
    self.resids_2 = resids[1]
    self.title = "Distance-difference matrix"

  def set_title (self, title) :
    self.title = title

def draw_plot (ddm,
               figure) :
  plot = figure.add_axes([0.1,0.1,0.8,0.8])
  im = plot.imshow(ddm.m.as_numpy_array(), origin="lower")
  plot.set_xlabel("Residue ID")
  plot.set_ylabel("Residue ID")
  xticklabels = []
  for x in plot.get_xticks() :
    if (x >= 0) and (x < len(ddm.resids_1)) :
      xticklabels.append(ddm.resids_1[int(x)])
    else :
      xticklabels.append("")
  plot.set_xticklabels(xticklabels)
  yticklabels = []
  for y in plot.get_yticks() :
    if (y >= 0) and (y < len(ddm.resids_2)) :
      yticklabels.append(ddm.resids_1[int(y)])
    else :
      yticklabels.append("")
  plot.set_yticklabels(yticklabels)
  cb = figure.colorbar(im, ax=plot, fraction=0.05, aspect=40)
  cb.set_label("Change in C-alpha:C-alpha distance")
  plot.set_title(ddm.title)

def display_plot_pylab (ddm, savefig=False) :
  from matplotlib import pyplot as plt
  figure = plt.figure(figsize=(12,10))
  draw_plot(ddm, figure)
  if (savefig) :
    figure.savefig("distance_difference.png", format="png")
  else :
    plt.show()

def validate_params (params) :
  params = params.calculate_matrix
  for i in [1,2] :
    file_param = getattr(params, "model_%d" % i)
    chain_param = getattr(params, "chain_%d" % i)
    if (file_param is None) :
      raise Sorry("Missing file for model_%d!" % i)
    elif (not os.path.isfile(file_param)) :
      raise Sorry("%s is not a file or does not exist." % file_param)
    if (chain_param is not None) and ((len(chain_param) > 2) or
        (len(chain_param) == 0)) :
      raise Sorry(("Invalid chain ID '%s' - must be 1 or 2 characters in "+
        "length if explicitly set.") % chain_param)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:], display_plot=True)
