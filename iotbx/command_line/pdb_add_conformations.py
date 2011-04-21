# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.add_conformations

import libtbx.phil
from libtbx import runtime_utils
from libtbx.utils import Sorry, Usage
import libtbx.load_env # import dependency
import string
import os
import sys

master_phil = libtbx.phil.parse("""
add_conformations
  .caption = This utility will duplicate any set of atoms (by default, the \
    entire input model) to create alternate conformations.  If the new \
    occupancy is not specified, it will be split evently between each \
    conformation.  Please be aware that if alternate conformations or \
    reduced-occupancy atoms are already present in the starting model, the \
    program behavior is not well-defined, and it may fail.
  .style = auto_align caption_img:icons/custom/iotbx.pdb.add_conformations64.png
{
  pdb_file = None
    .type = path
    .short_caption = PDB file
    .style = bold file_type:pdb noauto
  atom_selection = None
    .type = str
    .input_size = 400
    .style = bold noauto
  output = None
    .type = path
    .short_caption = Output file
    .style = bold file_type:pdb new_file noauto
  n_confs = 2
    .type = int
    .short_caption = Total number of conformations
    .style = bold spinner min:2 max:16 noauto
  new_occ = None
    .type = float
    .short_caption = New occupancy
    .style = noauto
  new_altloc = B
    .type = str
    .short_caption = Start at altloc
    .input_size = 64
    .style = noauto
}
""")
master_params = master_phil # XXX backwards compatibility for phenix gui

def run (args=(), params=None, out=sys.stdout) :
  if (len(args) == 0) and (params is None) :
    raise Usage("iotbx.pdb.add_conformations model.pdb [selection=...]\n"+
      "Full parameters:\n" + master_phil.as_str())
  from iotbx import file_reader
  pdb_in = None
  if (params is None) :
    user_phil = []
    interpreter = master_phil.command_line_argument_interpreter(
      home_scope="")
    for arg in args :
      if os.path.isfile(arg) :
        f = file_reader.any_file(os.path.abspath(arg))
        if (f.file_type == "pdb") :
          pdb_in = f.file_object
          user_phil.append(libtbx.phil.parse(
            "add_conformations.pdb_file=\"%s\"" % f.file_name))
        elif (f.file_type == "phil") :
          user_phil.append(f.file_object)
        else :
          raise Sorry("Unknown file type '%s' (%s)" % (f.file_type, arg))
      else :
        try :
          arg_phil = interpreter.process(arg=arg)
        except RuntimeError, e :
          raise Sorry("Error parsing '%s': %s" % (arg, str(e)))
        else :
          user_phil.append(arg_phil)
    params = master_phil.fetch(sources=user_phil).extract()
  validate_params(params)
  params = params.add_conformations
  if (pdb_in is None) :
    f = file_reader.any_file(params.pdb_file)
    f.assert_file_type("pdb")
    pdb_in = f.file_object
  if (params.new_occ is None) :
    params.new_occ = 1.0 / params.n_confs
    print >> out, "Setting new occupancy to %.2f" % params.new_occ
  from scitbx.array_family import flex
  hierarchy = pdb_in.construct_hierarchy()
  all_atoms = hierarchy.atoms()
  all_atoms.reset_i_seq()
  n_atoms = all_atoms.size()
  if (params.atom_selection is not None) :
    cache = hierarchy.atom_selection_cache()
    selection = cache.selection(params.atom_selection).as_int()
    if (selection.count(1) == 0) :
      raise Sorry("Empty selection.")
  else :
    selection = flex.int(n_atoms, 1)
  for i_seq, atom in enumerate(all_atoms) :
    atom.tmp = selection[i_seq]
  for model in hierarchy.models() :
    for chain in model.chains() :
      n_confs = len(chain.conformers())
      for residue_group in chain.residue_groups() :
        i_ag = 0
        for atom_group in residue_group.atom_groups() :
          old_atoms = atom_group.atoms()
          flags = old_atoms.extract_tmp_as_size_t()
          if (flags.count(1) > 0) :
            if (atom_group.altloc != "") :
              if (n_confs > 0) :
                atoms_err = [ a.format_atom_record() for a in old_atoms ]
                raise Sorry("Atom group included in selection already has one "+
                  "or more alternate conformers:\n" + "\n".join(atoms_err))
              elif (atom_group.altloc == params.new_altloc) :
                atoms_err = [ a.format_atom_record() for a in old_atoms ]
                raise Sorry("Atom group included in selection has an altloc "+
                  "identical to the new_altloc parameter:\n" +
                  "\n".join(atoms_err))
            else :
              atom_group.altloc = "A"
            old_occ = old_atoms.extract_occ()
            new_altloc = params.new_altloc
            for n in range(params.n_confs - 1) :
              new_group = atom_group.detached_copy()
              new_group.altloc = new_altloc
              if (flags.count(0) > 0) :
                j_atom = 0
                k_atom = 0
                for atom in new_group.atoms() :
                  if (atom.tmp == 0) :
                    new_group.remove_atom(j_atom)
                  else :
                    atom.set_occ(params.new_occ)
                    old_atoms[k_atom].set_occ(old_occ[k_atom] - params.new_occ)
                    j_atom += 1
                  k_atom += 1
              else :
                for old_atom, new_atom in zip(old_atoms, new_group.atoms()) :
                  old_occ = old_atom.occ - params.new_occ
                  if (old_occ == 0) :
                    print >> out, "WARNING: zero-occupancy atom:"
                    print >> out, old_atom.format_atom_record()
                  elif (old_occ < 0) :
                    raise Sorry("Atom occupancy dropped below zero:\n" +
                      old_atom.format_atom_record() + "\nnew_occ may be set "+
                      "too high.")
                  old_atom.set_occ(old_atom.occ - params.new_occ)
                  new_atom.set_occ(params.new_occ)
              assert (new_group.atoms().size() == flags.count(1))
              residue_group.insert_atom_group(i_ag + 1, new_group)
              i_ag += 1
              new_altloc = increment_altloc(new_altloc)
          i_ag += 1
  n_atoms_new = hierarchy.atoms().size()
  hierarchy.atoms().reset_i_seq()
  hierarchy.atoms_reset_serial()
  if (params.output is None) :
    base_name = os.path.basename(params.pdb_file)
    params.output = os.path.splitext(base_name)[0] + "_split.pdb"
  f = open(params.output, "w")
  f.write("\n".join(pdb_in.crystallographic_section()) + "\n")
  f.write(hierarchy.as_pdb_string())
  f.close()
  print >> out, "Old model: %d atoms" % n_atoms
  print >> out, "Modified model: %d atoms" % n_atoms_new
  print >> out, "Wrote %s" % params.output
  return params.output

def increment_altloc (altloc) :
  if altloc.isupper() :
    letters = string.uppercase
  elif altloc.islower() :
    letters = string.lowercase
  elif altloc.isdigit() :
    letters = string.digits
  else :
    raise Sorry("altloc must be a letter or digit.")
  i = letters.index(altloc) + 1
  if (i == len(letters)) :
    raise RuntimeError("Uh-oh, out of altlocs.")
  return letters[i]

def validate_params (params) :
  params = params.add_conformations
  if (params.pdb_file is None) :
    raise Sorry("Please specify a PDB file!")
  if (params.new_altloc is None) or (len(params.new_altloc) != 1) :
    raise Sorry("new_altloc must be a single character (e.g. 'B')")
  if (params.n_confs < 2) :
    raise Sorry("Number of conformations must be at least 2!")
  if (params.new_occ is None) :
    params.new_occ = 1.0 / params.n_confs
    #print >> out, "Setting new occupancy to %.2f" % params.new_occ
  if (params.new_occ < 0) or (params.new_occ > 1) :
    raise Sorry("new_occ must be between 0 and 1.0")
  return True

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    return run(args=list(self.args), out=sys.stdout)

if (__name__ == "__main__") :
  run(sys.argv[1:])
