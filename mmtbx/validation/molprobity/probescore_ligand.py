from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import easy_run
from libtbx.utils import Sorry
import sys

class probescore_result():
  def __init__(self, model_piece, probe_sel, id_str, nuclear):
    self.id_str = id_str
    self.probe_lines = None
    self.probe_lines = self.run_probescore(model_piece, probe_sel, nuclear)
    self.result = self.parse_probescore_lines(self.probe_lines)

  def run_probescore(self, pdb_string, sel_str, nuclear):
    #run probescore on a subset of the model
    if nuclear:
      nuclear_flag = "-nuclear"
    else:
      nuclear_flag = ""
    probe_command = 'phenix.probe -c -both %s "%s" "not (%s)" -' % (nuclear_flag,sel_str,sel_str)
    probe_out = easy_run.fully_buffered(probe_command, stdin_lines=pdb_string)
    if (probe_out.return_code != 0):
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    probe_lines = probe_out.stdout_lines
    return probe_lines

  def print_as_raw(self, out=sys.stdout):
    #print complete raw probe output
    print("\n----------------------------------------")
    print("Selection: %s" % self.id_str)
    for line in self.probe_lines:
      print(line, file=out)
    print("Selection: %s" % self.id_str)
    print("----------------------------------------")

  def print_as_digest(self, out=sys.stdout):
    #print a human-readable digest of the results
    ligand = self.result["ligand"]
    nonligand = self.result["nonligand"]
    if int(ligand["max_dots"]) == 0:
      ligand_dots_percent = 0
    else:
      ligand_dots_percent = float(ligand["contact_dots"])/float(ligand["max_dots"])*100
    if float(ligand["max_area"]) == 0:
      ligand_area_percent = 0
    else:
      ligand_area_percent = float(ligand["contact_area"])/float(ligand["max_area"])*100
    total_probescore = float(ligand["total_score"])+float(nonligand["total_score"])
    print("\nSelection: %s" % self.id_str)
    print("Selection dots: %s / %s = %.2f%%" % (ligand["contact_dots"],ligand["max_dots"],ligand_dots_percent), file=out)
    print("Selection contact area: %s / %s A^2 = %.2f%%" % (ligand["contact_area"],ligand["max_area"],ligand_area_percent), file=out)
    print("                     score  normalized", file=out)
    print("selection Hbond:  %8s %11s" % (ligand["hb_score"],ligand["hb_score_norm"]), file=out)
    print("selection Overlap:%8s %11s" % (ligand["overlap_score"],ligand["overlap_score_norm"]), file=out)
    print("selection Vdw:    %8s %11s" % (ligand["vdw_score"],ligand["vdw_score_norm"]), file=out)
    print("sum of contacts:  %8s %11s" % (ligand["total_score"],ligand["total_score_norm"]), file=out)
    print("\nprobescore:  %13s" % (total_probescore), file=out)
    print("\n")

  def print_as_oneline(self, filename="", out=sys.stdout):
    #print digest content as a single machine-readable line
    ligand = self.result["ligand"]
    nonligand = self.result["nonligand"]
    if int(ligand["max_dots"]) == 0:
      ligand_dots_percent = 0
    else:
      ligand_dots_percent = float(ligand["contact_dots"])/float(ligand["max_dots"])*100
    if float(ligand["max_area"]) == 0:
      ligand_area_percent = 0
    else:
      ligand_area_percent = float(ligand["contact_area"])/float(ligand["max_area"])*100
    total_probescore = float(ligand["total_score"])+float(nonligand["total_score"])
    print(":".join([filename,
                    self.id_str,
                    ligand["atoms"],
                    ligand["contact_dots"],
                    ligand["max_dots"],
                    "%.2f" % ligand_dots_percent,
                    ligand["contact_area"],
                    ligand["max_area"],
                    "%.2f" % ligand_area_percent,
                    ligand["hb_score"],
                    ligand["hb_score_norm"],
                    ligand["overlap_score"],
                    ligand["overlap_score_norm"],
                    ligand["vdw_score"],
                    ligand["vdw_score_norm"],
                    ligand["total_score"],
                    ligand["total_score_norm"],
                    "%.2f" % total_probescore]))

  def parse_probescore_lines(self, probescore_lines):
    ligand_section = False
    other_section = False
    for line in probescore_lines:
      #output has a 1->2 section for contacts from the ligand
      # and a 2->1 section for contacts from the rest of the structure.
      if line.startswith("subgroup: 1->2"):
        ligand_section = True
        other_section = False
      elif line.startswith("subgroup: 2->1"):
        ligand_section = False
        other_section = True
      #-----
      if ligand_section:
        if line.startswith("atoms"):
          ligand_atoms_selected = line.split()[2]
        elif line.startswith("potential dots"):
          ligand_potential_dots = line.split()[2]
        elif line.startswith("potential area"):
          ligand_potential_area = line.split()[2]
        #-----
        elif line.startswith("     tot contact"):
          ligand_vdw_score = line.split()[4]
          ligand_vdw_score_normalized  = line.split()[5]
        elif line.startswith("     tot overlap"):
          ligand_o_score = line.split()[4]
          ligand_o_score_normalized  = line.split()[5]
        elif line.startswith("     tot  H-bond"):
          ligand_hb_score = line.split()[4]
          ligand_hb_score_normalized  = line.split()[5]
        elif line.startswith("       grand tot"):
          ligand_total_dots = line.split()[2]
          ligand_total_score = line.split()[4]
          ligand_total_score_normalized  = line.split()[5]
        elif line.startswith("contact"):
          ligand_contact_area = line.split()[3]
      if other_section:
        #contacts from protein/na/other to the ligand
        #atom selection  is variable, based on "within",
        #  so normalization and area aren't very meaningful
        if line.startswith("atoms"):
          other_atoms_selected = line.split()[2]
        elif line.startswith("potential dots"):
          other_potential_dots = line.split()[2]
        elif line.startswith("potential area"):
          other_potential_area = line.split()[2]
        #-----
        elif line.startswith("     tot contact"):
          other_vdw_score = line.split()[4]
          #other_vdw_score_normalized  = line.split()[5]
        elif line.startswith("     tot overlap"):
          other_o_score = line.split()[4]
          #other_o_score_normalized  = line.split()[5]
        elif line.startswith("     tot  H-bond"):
          other_hb_score = line.split()[4]
          #other_hb_score_normalized  = line.split()[5]
        elif line.startswith("       grand tot"):
          other_total_dots = line.split()[2]
          other_total_score = line.split()[4]
          #other_total_score_normalized  = line.split()[5]
        elif line.startswith("contact"):
          other_contact_area = line.split()[3]

    ligand = {"atoms":ligand_atoms_selected,
              "max_dots":ligand_potential_dots,
              "contact_dots":ligand_total_dots,
              "max_area":ligand_potential_area,
              "contact_area":ligand_contact_area,
              "vdw_score":ligand_vdw_score,
              "vdw_score_norm":ligand_vdw_score_normalized,
              "hb_score":ligand_hb_score,
              "hb_score_norm":ligand_hb_score_normalized,
              "overlap_score":ligand_o_score,
              "overlap_score_norm":ligand_o_score_normalized,
              "total_score":ligand_total_score,
              "total_score_norm":ligand_total_score_normalized
              }
    nonligand = {"atoms":other_atoms_selected,
              "max_dots":other_potential_dots,
              "contact_dots":other_total_dots,
              "max_area":other_potential_area,
              "contact_area":other_contact_area,
              "vdw_score":other_vdw_score,
              #"vdw_score_norm":other_vdw_score_normalized,
              "hb_score":other_hb_score,
              #"hb_score_norm":other_hb_score_normalized,
              "overlap_score":other_o_score,
              #"overlap_score_norm":other_o_score_normalized,
              "total_score":other_total_score,
              #"total_score_norm":other_total_score_normalized
              }
    return {"ligand":ligand, "nonligand":nonligand}

class probescore():
  def __init__(self, model, selection_string_list, has_h, nuclear=False, out=sys.stdout):
    self.atoms = model.get_atoms()
    self.all_bsel = flex.bool(self.atoms.size(), False)
    self.results = []
    if not has_h:
      model_piece = self.add_hydrogens(model, nuclear=nuclear)
      ###TODO: read the Reduce'd version back in and do the model_piece thing on it
      ###This will be good for kinemage printing, if I get there
    for sel_str in selection_string_list:
      probe_sel, id_str = self.convert_to_probe_sel(model, sel_str)
      #probe_sel is the selection string rendered into probe syntax "chainA  100|chainA  350"
      #id_str is human readable for printing "A ASP 100,A SPD 350"
      if has_h:
        model_piece = self.get_model_piece(model, sel_str)
      self.results.append(probescore_result(model_piece, probe_sel, id_str, nuclear))

  def print_as_raw(self, out=sys.stdout):
    #print complete raw probe output
    for result in self.results:
      result.print_as_raw(out=out)

  def print_as_digest(self, out=sys.stdout):
    #print a human-readable digest of the results
    for result in self.results:
      result.print_as_digest(out=out)

  def print_as_oneline(self, filename=None, out=sys.stdout):
    #print digest content as a single machinereadable line
    header = ["input_file",
              "selection",
              "atoms",
              "contact_dots",
              "max_dots",
              "dots_percent",
              "contact_area",
              "max_area",
              "area_percent",
              "hb_score",
              "hb_score_norm",
              "overlap_score",
              "overlap_score_norm",
              "vdw_score",
              "vdw_score_norm",
              "total_score",
              "total_score_norm",
              "probescore"]
    print(":".join(header),file=out)
    for result in self.results:
      result.print_as_oneline(filename=filename, out=out)

  def convert_to_probe_sel(self, model, sel_str):
    #convert a Phenix selection string to a probe one
    #make the selection in Phenix,
    #then make a probe string enumerating every atom in the selection
    #instead, find the involved residues
    #this is gross, but is guarantees the same selection
    isel = model.iselection(string=sel_str)
    self.all_bsel.set_selected(isel, True)
    probe_residues = []
    id_strs = []
    for atom in self.atoms.select(isel):
      atomname = atom.name
      atomname = atomname.replace(" ","_") #probe uses _ chars for spaces in atom names for clarity
      resid = atom.parent().id_str()
      #ag.id_str() format: " TYR A  30 "
      altloc = resid[0:1]
      resname = resid[1:4]
      chain = resid[4:6].strip()
      resseq = resid[6:10]
      icode = resid[10:11].strip()
      probe_res_str = "chain%s %s%s" % (chain,resseq,icode)
      id_str = "%s %s%s%s" % (chain,resname,resseq,icode)
      probe_residues.append(probe_res_str)
      id_strs.append(id_str)
    unique_residues = set(probe_residues)
    unique_id_strs = set(id_strs)
    probe_sel = "|".join(unique_residues)
    id_str = ",".join(unique_id_strs)
    return probe_sel, id_str

  def add_hydrogens(self,model,nuclear=False):
    from mmtbx.validation.clashscore import check_and_add_hydrogen
    hierarchy = model.get_hierarchy()
    reduce_str, check = check_and_add_hydrogen(
      pdb_hierarchy=hierarchy,
      file_name=None,
      nuclear=nuclear,
      keep_hydrogens=True,
      verbose=False,
      model_number=0,
      n_hydrogen_cut_off=0,
      time_limit=120,
      allow_multiple_models=True,
      crystal_symmetry=None,
      do_flips=False,
      log=None)
    #There should be a check here for whether Reduce completed
    return reduce_str

  def get_model_piece(self, model, sel_str):
    #use "within" to get the relevant atoms from the model
    #convert those atoms to PDB string for into to probe
    #using subset of atoms should make probe faster for big models
    wide_selection = "within(%i, %s)" % (5, sel_str)
    isel = model.iselection(string=wide_selection)
    self.all_bsel.set_selected(isel, True)
    atom_list = []
    for atom in self.atoms.select(isel):
      atom_list.append(atom.format_atom_record())
    if not atom_list:
      raise Sorry('no atoms selected for "%s", please check selection' % sel_str)
    pdb_string = "\n".join(atom_list)
    return pdb_string

  def run_probescore(self, pdb_string, sel_str):
    #run probescore on a subset of the model
    probe_command = 'phenix.probe -c -both "%s" "not (%s)" -' % (sel_str,sel_str)
    probe_out = easy_run.fully_buffered(probe_command, stdin_lines=pdb_string)
    if (probe_out.return_code != 0):
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    result = probe_out.stdout_lines
    return result

