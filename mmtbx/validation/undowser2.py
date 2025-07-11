from __future__ import absolute_import, division, print_function

import sys
from mmtbx.programs import probe2
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.validation import validation, atoms, residue, atom_info
from mmtbx.validation.clashscore2 import check_and_add_hydrogen, probe_clashscore_manager
from mmtbx.validation.clashscore2 import remove_models_except_index
from libtbx.utils import null_out
import json
import iotbx.cli_parser
import mmtbx
import os
import tempfile

class ud_water(residue):
  __ud_water_attr__ = [
    "contacts",
    "model_id",
    "src_atom_id"
    #"max_b_factor",
  ]
  __slots__ = residue.__slots__ + __ud_water_attr__

  @staticmethod
  def header():
    return "%-20s %-20s %7s %7s %7s  %-20s" % ("Water ID", "Clashes with", "Water B", "Contact B", "Clash severity", "Category")

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if s != "contacts" and hasattr(self, s) ]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    wc_list = []
    for wc in self.contacts:
      wc_slots_list = [s for s in wc.__slots__ if s != "atoms_info" and s != "src_atom_id"]
      wc_slots_as_dict = {s: getattr(wc, s) for s in wc_slots_list if s != 'xyz' and s != "atom_selection"}
      wc_list.append(wc_slots_as_dict)
    slots_as_dict["water_contacts"] = wc_list
    #print({**slots_as_dict, **atom0_slots_as_dict})
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_string(self):
    return "%-20s %-20s  %7.3f" % (self.atoms_info[0].id_str(),
      self.atoms_info[1].id_str(), abs(self.overlap))

  def as_table_row_phenix(self):
    rows = []
    for wc in self.contacts:
      rows.append([self.id_str(), wc.atoms_info[1].id_str(), wc.src_b, wc.trg_b, wc.mingap, wc.get_category()])
    return rows

  def get_contacts(self):
    return self.contacts

#-u -q -mc -het -once -NOVDWOUT %s %s' % (probe_command, condensed_flag, nuclear_flag) "ogt%d not water" "ogt%d" -' % (ogt, ogt)
#       phenix.probe -4H -quiet -noticks -nogroup -dotmaster -mc -het -once -wat2wat 'water' 'water'

class undowserlyze(validation):
  __slots__ = validation.__slots__ + [
    "fast",
    "condensed_probe",
    "probe_file",
    "probe_undowser_manager",
    "water_count",
    "results_by_model"
  ]
  program_description = "Analyze waters for model"
  gui_list_headers = ["Water ID", "Clashes with", "Water B", "Contact B", "Clash severity", "Category"]
  gui_formats = ["%s", "%s", ".2f", ".2f", ".3f", "%s"]
  wx_column_widths = [150, 150, 150, 150, 150, 300] #actually set in GUI's Molprobity/Core.py
  html_header = """<html>
<head>
  <title>Summary table of water clashes</title>
</head>

<body>

<hr>
This table lists all HOH "waters" in the structure that have steric clashes. HOH are classified into common categories based on the atom they clash with.
<br><br>
A clashing HOH is very unlikely be be a real water, unless the clashing atom position is incorrect. The following categories provide guidance for correcting false HOH.
<br><br>
<b>Clash with polar</b> - HOH that clashes with polar groups may actually be a coordinated ion.
<br>
<b>Clash with nonpolar</b> - HOH that clashes with nonpolar groups may be a missing or displaced atom&ast;. Or it may be the first atom of an unmodeled alternate.
<br>
<b>Clash with both polar and nonpolar</b> - HOH that clashes with both polar and non-polar groups is unlikely to be an ion. If clashes are severe, a displaced atom is likely. If clashes and map are weak, the HOH may be entirely removable.
<br>
<b>Clash with water</b> - HOH-HOH clashes may be real waters that need to be modeled as alternates of compatible occupancy. Or they may be in the density of a sidechain alternate or a larger ligand.
<br>
<b>Clash with altloc</b> - HOH clashes involving one or more alternate conformations may be resolved by renaming some of the alternates.
<br><br>
<b>High B-factor</b> - HOH with clashes and minimal support in the map should be removed from the model. This table does not report map data directly, but a high B-factor is a likely warning sign that an HOH is a poor fit to the map.
<br>
<b>Severe clash</b> - HOH with severe clash overlap but good map support is likely to be a position where an atom is displaced.
<br><br>
&ast;<i>Displaced atom</i> indicates that a structural atom has been moved from its proper place in the model and replaced by HOH. Displaced sidechains are common. Moved atoms may be restored by local rebuilding.
<br>
<i>Missing atoms</i> have been entirely replaced by HOH. Removed atoms may be restored by modeling alternate conformations (especially sidechains), modeling ligands, or continuing a macromolecular mainchain.
<br><br>
These categories are general suggestions. Check your electron density; trust your intuition and experience. Prisant 2020 Prot Sci 29:315 (<a href="https://doi.org/10.1002/pro.3786">https://doi.org/10.1002/pro.3786</a>) illustrates 10 examples of clashing HOH cases.
<br>
<hr>
<br>
"""

  html_table = """
<br><br>
<hr>
<br>
<table border=1 width='100%'>
<tr bgcolor='#9999cc'><td rowspan='1' align='center'>Water ID</td>
<td align='center'>Clashes with</td>
<td align='center'>Water B</td>
<td align='center'>Contact B</td>
<td align='center'>Clash<br>Severity</td>
<td align='center'>Clash with Polar<br><small>May be ion</small></td>
<td align='center'>Clash with non-polar<br><small>Unmodeled alt or noise</small></td>
<td align='center'>Clash with water<br><small>Occ &lt;1 or ligand</small></td>
<td align='center'>Clash with altloc<br><small>Add or rename alts</small></td></tr>
"""

  def get_result_class(self):
    return ud_water

  def __init__(self,
      probe_parameters,
      data_manager,
      keep_hydrogens=True,
      nuclear=False,
      outliers_only=False,
      force_unique_chain_ids=False,
      b_factor_cutoff=None,
      save_modified_hierarchy=False,
      verbose=False,
      do_flips=False,
      out=sys.stdout):
    validation.__init__(self)
    if verbose:
      if not nuclear:
        print("\nUsing electron cloud x-H distances and vdW radii")
      else:
        print("\nUsing nuclear cloud x-H distances and vdW radii")
    import iotbx.pdb
    from scitbx.array_family import flex
    from mmtbx.validation import utils

    data_manager_model = data_manager.get_model()
    # Fix up bogus unit cell when it occurs by checking crystal symmetry.
    # @todo reduce_hydrogens.py:run() says: TODO temporary fix until the code is moved to model class
    cs = data_manager_model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      data_manager_model = shift_and_box_model(model = data_manager_model)

    # If we've been asked to, add hydrogens to all of the models in the PDB hierarchy
    # associated with our data_manager_model.
    data_manager_model,_ = check_and_add_hydrogen(
      probe_parameters=probe_parameters,
      data_manager_model=data_manager_model,
      nuclear=nuclear,
      verbose=verbose,
      keep_hydrogens=keep_hydrogens,
      do_flips = do_flips,
      log=null_out())

    # First we must rebuild the model from the new hierarchy so that the copy can succeed.
    # Make a copy of the original model to use for submodel processing, we'll trim atoms out
    # of it for each submodel.
    data_manager_model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = data_manager_model.get_hierarchy(),
      stop_for_unknowns = False,
      crystal_symmetry  = data_manager_model.crystal_symmetry(),
      restraint_objects = None,
      log               = null_out())
    original_model = data_manager_model.deep_copy()

    pdb_hierarchy = data_manager_model.get_hierarchy()
    n_models = len(pdb_hierarchy.models())
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=pdb_hierarchy)

    # Get information about the waters in the model
    sel_cache = pdb_hierarchy.atom_selection_cache()
    water_sel= sel_cache.selection("water")
    hierarchy_waters = pdb_hierarchy.select(water_sel)
    self.water_count = len(list(hierarchy_waters.residue_groups()))
    water_xyzs = {}
    self.results_by_model = {}
    for atom in hierarchy_waters.atoms():
      if atom.name == " O  ":
        water_xyzs[atom.parent().parent().resid()] = atom.xyz
        model_id = atom.parent().parent().parent().parent().id
        if model_id not in self.n_total_by_model:
          self.n_total_by_model[model_id] = 0
        self.n_total_by_model[model_id] += 1

    for i_mod, model in enumerate(pdb_hierarchy.models()):
      if model.id not in self.results_by_model:
        self.results_by_model[model.id] = []
        self.n_outliers_by_model[model.id] = 0
      # Select only the current submodel from the hierarchy
      submodel = original_model.deep_copy()
      remove_models_except_index(submodel, i_mod)

      # Construct a hierarchy for the current submodel
      r = iotbx.pdb.hierarchy.root()
      mdc = submodel.get_hierarchy().models()[0].detached_copy()
      r.append_model(mdc)

      occ_max = flex.max(r.atoms().extract_occ())

      # Make yet another model for the new hierarchy
      subset_model_manager = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = r,
        stop_for_unknowns = False,
        crystal_symmetry  = submodel.crystal_symmetry(),
        restraint_objects = None,
        log               = null_out())

      self.probe_undowser_manager = probe_undowser_manager(
        nuclear=nuclear,
        verbose=verbose,
        model_id=model.id)
      self.probe_undowser_manager.run_probe_undowser(data_manager, subset_model_manager)
      water_contacts = self.probe_undowser_manager.get_water_contacts()
      for src_atom_id, wc_list in water_contacts.items():
        src_chain_id=probe_undowser_manager.get_chain_probe_atom_id(src_atom_id)
        src_resseq=probe_undowser_manager.get_resseq_probe_atom_id(src_atom_id)
        src_icode=probe_undowser_manager.get_icode_probe_atom_id(src_atom_id)
        src_resname=probe_undowser_manager.get_resname_probe_atom_id(src_atom_id)
        src_altloc=probe_undowser_manager.get_altloc_probe_atom_id(src_atom_id)
        #selection_string = "model {} and chain {} and resseq '{}' and altid '{}' and name O".format(model.id, src_chain_id, src_resseq, src_altloc)
        result = ud_water(
          outlier=True,
          src_atom_id=src_atom_id,
          model_id=model.id,
          chain_id=src_chain_id,
          resseq=src_resseq,
          icode=src_icode,
          resname=src_resname,
          altloc=src_altloc,
          segid=None, # XXX ???
          xyz=water_xyzs[src_resseq+src_icode],
          contacts=wc_list,
        )
        #if (not outliers_only or is_outlier):
        self.results.append(result)
        self.results_by_model[model.id].append(result)
        self.n_outliers_by_model[model.id] += 1

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "undowser"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results

    for model_id in self.n_total_by_model.keys():
      summary_results[model_id] = { "num_outliers" : self.n_outliers_by_model[model_id],
        "num_waters" : self.n_total_by_model[model_id],
      }
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def as_HTML(self):
    return self.as_html_table()

  #defaults to returning table corresponding to the first model, even if model_id isn't a valid input
  def as_html_table(self, model_id=""):
    if model_id not in self.results_by_model:
      model_id = sorted(self.results_by_model.keys())[0]
    model_results = self.results_by_model[model_id]
    #this could probably be more efficient, this recreates a "water_contacts" object to avoid changing CJW's original HTML code
    water_contacts = {}
    for ud_water in model_results:
      water_contacts[ud_water.src_atom_id]=ud_water.get_contacts()

    #water_contacts = self.probe_undowser_manager.get_water_contacts()
    html_string = self.html_header
    if self.water_count == 0:
      html_string = html_string+"SUMMARY: %i waters out of %i have clashes (%.2f%%)" % (0, 0, 0)
    else:
      html_string = html_string+"SUMMARY: %i waters out of %i have clashes (%.2f%%)" % (len(water_contacts), self.n_total_by_model[model_id], len(water_contacts)/self.n_total_by_model[model_id]*100)
    html_string = html_string + self.html_table

    contact_keys = sorted(water_contacts, key=lambda c: (-1*cumulative_severity(c, water_contacts)))
    #simple reverse sorting may put waters with the same cumulative severity in reverse sequence order

    row_number = 0
    row_color = ['#eaeaea','#ffffff']
    for contact_key in contact_keys:
      water = water_contacts[contact_key]
      water.sort(key=lambda w: (-1*float(w.mingap)))
      bgcolor = row_color[row_number%2]
      html_string = html_string+"<tr bgcolor=%s><td rowspan='%i' ><pre><code>%s</code></pre></td>\n" % (bgcolor,len(water),water[0].format_src_id_str())
      row_number+=1
      clashcount = 0
      for clash in water:
        if clashcount: html_string = html_string+'<tr bgcolor=%s>' % bgcolor
        clashcount+=1
        html_string = html_string+"<td><pre><code>%s</code></pre></td>%s%s<td bgcolor='%s'>%s</td>%s%s%s%s</tr>\n" % (clash.format_trg_id_str(), clash.write_b_cell(clash.src_b), clash.write_b_cell(clash.trg_b), clash.color_clash_severity(), clash.mingap, clash.write_polar_cell(), clash.write_nonpolar_cell(), clash.write_other_water_cell(), clash.write_altloc_cell())
    html_string = html_string+"</table>"
    return html_string
      #out.write("<td>%s</td><td bgcolor='%s'>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" % (clash.format_trg_id_str(), clash.color_clash_severity(), clash.mingap,  clash.does_it_clash_with_altloc(), clash.does_it_clash_with_polar(), clash.does_it_clash_with_nonpolar(), clash.does_it_clash_with_other_water()))
      #out.write("<td>%s</td><td bgcolor='%s'>%s</td><td %s></td><td %s></td><td %s></td><td %s></td></tr>\n" % (clash.format_trg_id_str(), clash.color_clash_severity(), clash.mingap, color_cell(clash.does_it_clash_with_polar()), color_cell(clash.does_it_clash_with_nonpolar()), color_cell(clash.does_it_clash_with_other_water()), color_cell(clash.does_it_clash_with_altloc()) ))

    #tr-level bgcolors for alternating table stripes
    #bgcolor='#9999cc' (blue for column headers)
    #bgcolor='#ffffff' (white background)
    #bgcolor='#f0f0f0' (original gray for stripes)
    #bgcolor='#eaeaea' (slightly darker gray for stripes)

class water_contact(atoms):
  wc_slots = ["model_id",
    "src_atom_id",
    "src_chain",
    "src_resseq",
    "src_icode",
    "src_resname",
    "src_atom",
    "src_altloc",
    "trg_heavy_atom",
    "trg_chain",
    "trg_resseq",
    "trg_icode",
    "trg_resname",
    "trg_atom",
    "trg_altloc",
    "mingap",
    "src_b",
    "trg_b",
    "category",
    ]

  __slots__ = atoms.__slots__ + wc_slots

  def __init__(self, model_id, src_atom_id, trg_atom_id, trg_heavy_atom, mingap, src_b, trg_b):
    self.model_id = model_id
    self.src_atom_id = src_atom_id
    #self.trg_atom_id = trg_atom_id
    self.trg_heavy_atom = trg_heavy_atom
    # A2778 HOH  O  A
    self.atoms_info = []
    #ccnnnnirrr?aaaal
    self.src_chain = probe_undowser_manager.get_chain_probe_atom_id(src_atom_id)
    self.src_resseq = probe_undowser_manager.get_resseq_probe_atom_id(src_atom_id)
    self.src_icode = probe_undowser_manager.get_icode_probe_atom_id(src_atom_id)
    self.src_resname = probe_undowser_manager.get_resname_probe_atom_id(src_atom_id)
    self.src_atom = probe_undowser_manager.get_atom_probe_atom_id(src_atom_id)
    self.src_altloc = probe_undowser_manager.get_altloc_probe_atom_id(src_atom_id)
    self.atoms_info.append(atom_info(
      model_id=model_id,
      chain_id=self.src_chain,
      resseq=self.src_resseq,
      icode=self.src_icode,
      resname=self.src_resname,
      altloc=self.src_altloc,
      name=self.src_atom)
    )

    self.trg_chain = probe_undowser_manager.get_chain_probe_atom_id(trg_atom_id)
    self.trg_resseq = probe_undowser_manager.get_resseq_probe_atom_id(trg_atom_id)
    self.trg_icode = probe_undowser_manager.get_icode_probe_atom_id(trg_atom_id)
    self.trg_resname = probe_undowser_manager.get_resname_probe_atom_id(trg_atom_id)
    self.trg_atom = probe_undowser_manager.get_atom_probe_atom_id(trg_atom_id)
    self.trg_altloc = probe_undowser_manager.get_altloc_probe_atom_id(trg_atom_id)
    self.atoms_info.append(atom_info(
      model_id=model_id,
      chain_id=self.trg_chain,
      resseq=self.trg_resseq,
      icode=self.trg_icode,
      resname=self.trg_resname,
      altloc=self.trg_altloc,
      name=self.trg_atom)
    )

    self.mingap = mingap.lstrip('-')
    #main table lists all overlaps as positive
    #the simple strip works as long as we only look at clashes

    self.src_b = float(src_b)
    self.trg_b = float(trg_b)
    self.outlier = True
    self.score = 0.0
    self.category = self.get_category()

  def get_category(self):
    category = "uncategorized"
    if self.does_it_clash_with_polar():
      category = "polar clash"
    if self.does_it_clash_with_nonpolar():
      category = "nonpolar clash"
    if self.does_it_clash_with_other_water():
      category = "water clash"
    if self.does_it_clash_with_altloc():
      if category != "uncategorized":
        category = category + ", altloc clash"
      else:
        category = "altloc clash"
    return category

  def trg_is_H(self):
    if self.trg_atom.strip().startswith('H'): return True
    else: return False

  def trg_is_charged_N(self):
  #what N's can be charged in nucleic acids?
  #how to handle the great variety of ligands?
    if self.trg_resname == 'LYS' and self.trg_atom == ' NZ ':
      return True
    elif self.trg_resname == 'ARG' and self.trg_atom in [' NH1',' NH2']:
      return True
    elif self.trg_resname == 'HIS' and self.trg_atom in [' ND1',' NE2']:
      return True
    else:
      return False

  def format_src_id_str(self):
    return ":".join([self.src_chain,self.src_resseq+self.src_icode,self.src_resname,self.src_altloc])

  def format_trg_id_str(self):
    #return self.trg_heavy_atom+" of "+":".join([self.trg_chain,self.trg_resseq+self.trg_icode,self.trg_resname,self.trg_altloc])
    return self.trg_atom+" of "+":".join([self.trg_chain,self.trg_resseq+self.trg_icode,self.trg_resname,self.trg_altloc])

  def color_clash_severity(self):
    mingap = float(self.mingap)
    if mingap < 0.5: return '#ffb3cc'
    elif mingap > 0.9: return '#ee4d4d'
    else: return '#ff76a9'

  def does_it_clash_with_polar(self):
    #Polar atoms in proteins are O, N, and S
    #Nucleic acids also have P, but that appears completely shielded by O's
    #All O are polar.  Are any N sufficiently non-polar that they wouldn't coordinate with metal?
    #Another water does not count as polar for these purposes
    #Not currently considering what polar atoms might be in ligands
    if self.trg_resname == "HOH":
      return ''
    if self.trg_heavy_atom in ['O','S']:
      return True
    if self.trg_heavy_atom == 'N':
      if self.trg_is_H():
        return True #H on N is polar
      elif self.trg_is_charged_N():
        return True
    return ''

  def does_it_clash_with_altloc(self):
    #Goal is to find any clashes that might be resolved by different altloc naming
    #A-A clashes, A-_ clashes, and _-A clashes
    #Probe should automatically ignore any A-B contacts, since those are in fully different alts
    if self.src_altloc.strip() or self.trg_altloc.strip():
      return True
    else:
      return ''

  def does_it_clash_with_other_water(self):
    if self.trg_resname == "HOH":
      return True
    else:
      return ''

  def does_it_clash_with_nonpolar(self):
    #This is the most complex one and will likely receive iterations
    #Start with simple check for non-polars
    if self.trg_heavy_atom in ['C']:
      return True
    elif self.trg_is_H():
      if self.trg_heavy_atom in ['O','N','S']:
        return ''
      else:
        return True
    elif self.trg_heavy_atom in ['N']:
      if not self.trg_is_charged_N():
        return True
    return ''

  def write_polar_cell(self):
    if self.does_it_clash_with_polar():
      trg_element = self.trg_atom.strip()[0:1]
      if trg_element == 'H':
        return "<td align='center' bgcolor='%s'>&minus; ion</td>" % self.color_clash_severity()
      elif trg_element in ['O','S']:
        return "<td align='center' bgcolor='%s'>&plus; ion</td>" % self.color_clash_severity()
      elif trg_element == 'N':
        return "<td align='center' bgcolor='%s'>&minus; ion</td>" % self.color_clash_severity()
      else:
        return "<td align='center' bgcolor='%s'>*</td>" % self.color_clash_severity()
    else:
      return "<td></td>"

  def write_nonpolar_cell(self):
    if self.does_it_clash_with_nonpolar():
      return "<td align='center' bgcolor='%s'>&times;</td>" % self.color_clash_severity()
    else:
      return "<td></td>"

  def write_altloc_cell(self):
    if self.does_it_clash_with_altloc():
      if self.does_it_clash_with_other_water():
        return "<td align='center' bgcolor='%s'>alt water</td>" % self.color_clash_severity()
      elif self.src_altloc.strip() and self.trg_altloc.strip():
        return "<td align='center' bgcolor='%s'>alt both sides</td>" % self.color_clash_severity()
      elif self.src_altloc.strip():
        return "<td align='center' bgcolor='%s'>alt water</td>" % self.color_clash_severity()
      else:
        return "<td align='center' bgcolor='%s'>alt partner atom</td>" % self.color_clash_severity()
      #return "<td align='center' bgcolor='%s'>&times;</td>" % self.color_clash_severity()
    else:
      return "<td></td>"

  def write_other_water_cell(self):
    if self.does_it_clash_with_other_water():
      return "<td align='center' bgcolor='%s'>&times;</td>" % self.color_clash_severity()
    else:
      return "<td></td>"

  def write_b_cell(self, b):
    return "<td>%.2f</td>" % b

def cumulative_severity(contact_key, water_contacts):
  cumulative_severity = 0
  water = water_contacts[contact_key]
  for clash in water:
    weighted_severity = float(clash.mingap) - 0.2
    cumulative_severity += weighted_severity
  return cumulative_severity

def color_cell(clash_check):
  if clash_check: return 'bgcolor=#ff76a9'
  else: return ''

#this subclasses probe_clashscore_manager but doesn't really use most of it. It's mainly to reuse the code for checking for probe existance
class probe_undowser_manager(probe_clashscore_manager):
  def __init__(self,
               nuclear=False,
               verbose=False,
               model_id=""):
    super().__init__(
               fast = False,
               condensed_probe=False,
               nuclear=nuclear,
               largest_occupancy=10,
               b_factor_cutoff=None,
               use_segids=False,
               verbose=verbose,
               model_id=model_id)
    self.water_contacts = {}

  def run_probe_undowser(self, data_manager, hydrogenated_model):

    # Construct override parameters and then run probe2 using them and delete the resulting
    # temporary file.
    tempName = tempfile.mktemp()
    parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
    #self.probe_command = "%s -u -q -mc -het -con -once -wat2wat -stdbonds -onlybadout 'water' 'all' -" % (self.probe_command)
    args = [
      "source_selection='water'".format(self.occupancy_frac),
      "target_selection='all'".format(self.occupancy_frac),
      "use_neutron_distances={}".format(self.nuclear),
      "approach=once",
      "include_water_water=True",
      "output.filename='{}'".format(tempName),
      "output.format=raw",
      "output.condensed=True",
      "output.report_vdws=False",
      "output.report_hydrogen_bonds=False",
      "output.only_report_bad_clashes=True",
      "ignore_lack_of_explicit_hydrogens=True",
    ]
    parser.parse_args(args)
    p2 = probe2.Program(data_manager, parser.working_phil.extract(),
                       master_phil=parser.master_phil, logger=null_out())
    p2.overrideModel(hydrogenated_model)
    dots, output = p2.run()
    probe_unformatted = output.splitlines()
    os.unlink(tempName)
    #print("probe unformatted:"+str(probe_unformatted))

    for line in probe_unformatted:
#>>name:pat:type:srcAtom:targAtom:dot-count:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
#SelfIntersect
#:1->2:bo: A2778 HOH  O  A: A 464 ARG  HD3 : 8:-0.581:-0.475: 36.509: 0.601: 18.650:0.238:-0.1485:O:C:36.622:0.786:18.552:30.97:17.73
#:1->2:bo: A2001 HOH  O  A: A2013 HOH  O  A:14:-0.792:-0.442:-12.858:17.914:-23.935:0.221:-0.1381:O:O:-12.726:18.090:-23.907:21.07:14.91
      n = line.split(':')
      contact_type = n[2]
      src_atom_id = n[3]
      trg_atom_id = n[4]
      dotcount = n[5]
      mingap = n[6]
      #src_heavy_atom = n[13] #always the water O
      trg_heavy_atom = n[14] #parent heavy atom of a clashing H
      #gap = n[7] #meaningless in -condensed mode, use mingap instead
      src_b = n[18]
      trg_b = n[19].strip()
      clash = water_contact(self.model_id, src_atom_id, trg_atom_id, trg_heavy_atom, mingap, src_b, trg_b)

      #sys.stdout.write(line)
      #sys.stdout.write('\n')
      polar_clash = clash.does_it_clash_with_polar()
      alt_clash = clash.does_it_clash_with_altloc()
      water_clash = clash.does_it_clash_with_other_water()
      nonpolar_clash = clash.does_it_clash_with_nonpolar()

      if src_atom_id not in self.water_contacts:
        self.water_contacts[src_atom_id] = []
      self.water_contacts[src_atom_id].append(clash)

  @staticmethod
  def get_chain_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[0:2].strip()

  @staticmethod
  def get_resseq_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[2:6]

  @staticmethod
  def get_icode_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[6:7]

  @staticmethod
  def get_resname_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[7:10]

  @staticmethod
  def get_atom_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[11:15]

  @staticmethod
  def get_altloc_probe_atom_id(atom_id):
    if len(atom_id) == 16:
      return atom_id[15:16]

  def get_water_contacts(self):
    return self.water_contacts
