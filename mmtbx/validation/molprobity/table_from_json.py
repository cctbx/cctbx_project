from __future__ import division
import os, sys
import json

red = '#ff9999'
yellow = '#ffff99'
green = '#99ff99'
MASTER_ORDER = ["clashscore",
                "ramalyze",
                "rotalyze",
                "cbetadev",
                "cablam",
                "rna_puckers",
                "rna_suites",
                "mp_bonds",
                "mp_angles",
                "omegalyze"]

def make_reskey(model_id, chain_id, resseq, icode):
  #create a unique string residue key with mmCIF-like placeholders for blank fields
  #The residue key is space-delimited
  return ' '.join([model_id.strip() or '.', chain_id.strip() or '.', resseq.strip(), icode.strip() or '?'])

class merged_validations():
  def __init__(self, hierarchy):
    self.data = {}
    self.indices = [] #reskeys in the same order as hierarchy
    self.validation_types = []
    self.summaries = {}
    self.model_ids = []
    for model in hierarchy.models():
      model_id = model.id.strip()
      self.model_ids.append(model_id)
      for chain in model.chains():
        chain_id = chain.id.strip()
        for rg in chain.residue_groups():
          #alts with different resnames?
          resseq = rg.resseq.strip()
          icode = rg.icode.strip()
          reskey = make_reskey(model_id, chain_id, resseq, icode)
          self.indices.append(reskey)
          self.data[reskey] = residue_bootstrap(rg, model_id, chain_id, resseq, icode)

  def add_summaries(self, json_list):
    for val_json in json_list:
      self.add_summary(val_json)

  def add_summary(self, val_json):
    val_type = val_json["validation_type"]
    self.summaries[val_type] = val_json["summary_results"]

  def add_validations(self, json_list):
    for val_json in json_list:
      val_type = val_json["validation_type"]
      #some validations especially those with multiple results per residue require special loading
      #these are clashscore, bond lengths, bond angles, chirals
      self.validation_types.append(val_type)
      if val_type == "clashscore":
        self.add_clashscore(val_json)
      elif val_type in ["mp_bonds", "mp_angles"]:
        self.add_geometry(val_json)
      else:
        self.add_validation(val_json)
    for reskey in self.data:
      #print(reskey)
      self.data[reskey].get_alternates()

  def add_validation(self, val_json):
    val_type = val_json["validation_type"]
    for result in val_json["flat_results"]:
      reskey = make_reskey(result["model_id"], result["chain_id"], result["resseq"].strip(), result["icode"])
      if not self.data[reskey].validations.get(val_type):
        self.data[reskey].validations[val_type] = {}
      self.data[reskey].validations[val_type][result["altloc"].strip()] = result
      if result["outlier"]:
        self.data[reskey].has_outlier = True

  def add_clashscore(self, val_json):
    assert val_json["validation_type"] == "clashscore"
    #each contact only gets one clash entry, so it needs to be stired twice for the residue table
    #  once for the source and once for the target
    #Usage of clash info should check whether reskey == srcreskey
    #  if not, then treat the target as the source
    for result in val_json["flat_results"]:
      srcreskey = make_reskey(result["model_id"], result["chain_id"], result["resseq"], result["icode"])
      if not self.data[srcreskey].validations.get("clashscore"):
        self.data[srcreskey].validations["clashscore"] = {}
      if not self.data[srcreskey].validations["clashscore"].get(result["altloc"].strip()):
        self.data[srcreskey].validations["clashscore"][result["altloc"].strip()] = []
      self.data[srcreskey].validations["clashscore"][result["altloc"].strip()].append(result)
      self.data[srcreskey].has_outlier = True

      trgreskey = make_reskey(result["model_id"], result["target_chain_id"], result["target_resseq"], result["target_icode"])
      if not self.data[trgreskey].validations.get("clashscore"):
        self.data[trgreskey].validations["clashscore"] = {}
      if not self.data[trgreskey].validations["clashscore"].get(result["target_altloc"].strip()):
        self.data[trgreskey].validations["clashscore"][result["target_altloc"].strip()] = []
      self.data[trgreskey].validations["clashscore"][result["target_altloc"].strip()].append(result)
      self.data[trgreskey].has_outlier = True
    #once all results are loaded, sort from worst to mildest clash
    for reskey in self.data:
      if not self.data[reskey].validations.get("clashscore"):
        continue
      for altloc in self.data[reskey].validations["clashscore"]:
        #clash overlaps are stored as negative values (vdw gaps would be positive)
        #so default ascending sort should put worst clash first
        self.data[reskey].validations["clashscore"][altloc].sort(key=lambda x: x["overlap"])

  def add_geometry(self, val_json):
    val_type = val_json["validation_type"]
    assert val_type in ["mp_bonds", "mp_angles"]
    for result in val_json["flat_results"]:
      #for both bonds and angles, the measure is assigned to the residue of the second atom (index 1) in the list
      #this only really matters at residue-residue joins
      reskey = make_reskey(result["atoms_model_id"][1], result["atoms_chain_id"][1], result["atoms_resseq"][1], result["atoms_icode"][1])

      #a bond or angle should be a mix of '' and possible a single alternate id
      #if multiple alternate ids are involved in a single bond, this will return the "highest" alternate
      #That's a weird case, and whether it shows up is based on how the hierarchy processes complex alternates.
      altloc = sorted(result["atoms_altloc"], reverse=True)[0].strip()

      if not self.data[reskey].validations.get(val_type):
        self.data[reskey].validations[val_type] = {}
      if not self.data[reskey].validations[val_type].get(altloc):
        self.data[reskey].validations[val_type][altloc] = []
      self.data[reskey].validations[val_type][altloc].append(result)
      self.data[reskey].has_outlier = True #needs to be outliers_only input
    #once all results are loaded, sort from worst to mildest
    #deviations can be too large or too small, so abs() is necessary
    for reskey in self.data:
      if not self.data[reskey].validations.get(val_type):
        continue
      for altloc in self.data[reskey].validations[val_type]:
        self.data[reskey].validations[val_type][altloc].sort(key=lambda x: abs(x["score"]), reverse=True)

  def get_table_order(self):
    #jsons may be received in any order
    table_order = []
    for x in MASTER_ORDER:
      if x in self.validation_types:
        table_order.append(x)
    return table_order

  def make_summary_table(self, model, red = '#ff9999', yellow = '#ffff99', green = '#99ff99'):
  # TODO: add a cleverness to do the first model if model==None
    #red = '#ff9999'
    #yellow = '#ffff99'
    #green = '#99ff99'
    lines = []
    table_order = self.get_table_order()
    if model.strip():
      lines.append("<br><b>Summary statistics: Model {model}</b>".format(model))
    else:
      lines.append("<br><b>Summary statistics</b>")
    lines.append("<table width='100%' cellspacing='1' border='2'>")
    lines += self.summary_table_clashscore(model, red, yellow, green)
    lines += self.summary_table_proteins(model, red, yellow, green)
    if "omegalyze" in table_order:
      lines = lines + self.summary_table_peptide_omegas(model, show_all=True)
    lines = lines + self.summary_table_nucleic_acids(model)
    lines = lines + self.summary_table_low_res(model)
    lines = lines + self.summary_table_additional(model)
    lines.append("</table>")
    lines.append(self.under_table_text(table_order))
    return "\n".join(lines)

  def under_table_text(self, table_order):
    notes = []
    notes.append("<small>In the two column results, the left column gives the raw count, right column gives the percentage.</small>")
    notes.append("<br><small>Some validations count different numbers of residues. E.g. Ramachandran's &phi; or &psi; is incomplete at chain ends, and Gly residues have no sidechains for rotamers.</small>")
    if "clashscore" in table_order:
      notes.append("<br><small>100<sup>th</sup> percentile is the best among structures of comparable resolution; 0<sup>th</sup> percentile is the worst.  For clashscore the comparative set of structures was selected in 2004, for MolProbity score in 2006.</small>")
      if "ramalyze" in table_order and "rotalyze" in table_order:
        notes.append("<br><small>MolProbity score combines the clashscore, rotamer, and Ramachandran evaluations into a single score, normalized to be on the same scale as X-ray resolution.</small>")
    if "rna_suites" in table_order:
      notes.append("<br><small>RNA backbone sugar-to-sugar suites are rotameric.  Outliers are RNA suites that don't fall into recognized rotamers.</small>")
    return "\n".join(notes)

  def summary_table_clashscore(self, model, red, yellow, green, resolution=None):
    resolution=2.0
    lines = []
    if "clashscore" in self.summaries:
      from mmtbx.validation.molprobity import percentile_clashscore
      data = self.summaries["clashscore"][model]
      pct_stats = percentile_clashscore.get_percentile_for_clashscore(data["clashscore"], resolution=resolution)
      #pct_stats = {"minres":minres, "maxres":maxres, "count":nSamples, "percentile":pctRank}
      if pct_stats["percentile"] >= 66: color=green
      elif pct_stats["percentile"] < 33: color=red
      else: color=yellow
      #ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
      if not resolution:
        range_text = "(N={count}, all resolutions)".format(count=pct_stats["count"])
      else:
        range_text = "(N={count}, {minres}&Aring; - {maxres}&Aring;)".format(count=pct_stats["count"], minres=pct_stats["minres"], maxres=pct_stats["maxres"])
      lines.append("""  <tr>
    <td align='center'>All-Atom Contacts</td>
    <td>Clashscore, all atoms</td>
    <td colspan='2' bgcolor='{color}'>&nbsp;{clashscore:.2f}</td>
    <td>{percentile}<sup>th</sup> percentile {range_text}</td>
  </tr>""".format(color=color, clashscore=data["clashscore"], percentile=pct_stats["percentile"], range_text=range_text))
    return lines

  def summary_table_proteins(self, model, red, yellow, green, resolution=None):
    protein_lines = []
    rows = 0
    if "rotalyze" in self.summaries:
      rows +=2
      data = self.summaries["rotalyze"][model]
      if data["outlier_percentage"] <= 0.3: color=green
      elif data["outlier_percentage"] > 1.5: color=red
      else: color=yellow
      protein_lines.append("  <tr>")
      protein_lines.append("""    <td>Poor rotamers</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{outlier_percentage:.2f}%</td>
    <td>Goal: &lt;0.3%</td>
  </tr>""".format(color= color, num_outliers=data["num_outliers"], num_residues=data["num_residues"], outlier_percentage=data["outlier_percentage"]))

      if data["num_residues"] == 0:
        color=green #none would be more accurate; does 'transparent' work?
        favored_pct = 0
      else:
        favored_pct = data["num_favored"]/data["num_residues"]*100.0
        if favored_pct >= 98: color=green
        elif favored_pct < 95: color=red
        else: color=yellow
      protein_lines.append("""  <tr>
    <td>Favored rotamers</td>
    <td bgcolor='{color}'>&nbsp;{num_favored} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{favored_pct:.2f}%</td>
    <td>Goal: &gt;98%</td>
  </tr>""".format(color=color, num_favored=data["num_favored"], num_residues=data["num_residues"], favored_pct=favored_pct))

    if "ramalyze" in self.summaries:
      rows += 2
      data = self.summaries["ramalyze"][model]
      if data["outlier_percentage"] <= 0.05: color=green
      elif data["outlier_percentage"] > 0.5 and data["num_outliers"] >= 2: color=red
      else: color=yellow
      protein_lines.append("  <tr>")
      protein_lines.append("""    <td>Ramachandran outliers</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{outlier_percentage:.2f}%</td>
    <td>Goal: &lt;0.05%</td>
  </tr>""".format(color=color, num_outliers=data["num_outliers"], num_residues=data["num_residues"], outlier_percentage=data["outlier_percentage"]))

      if data["favored_percentage"] >= 98: color=green
      elif data["favored_percentage"] < 95: color=red
      else: color=yellow
      protein_lines.append("""  <tr>
    <td>Ramachandran favored</td>
    <td bgcolor='{color}'>&nbsp;{num_favored} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{favored_percentage:.2f}%</td>
    <td>Goal: &gt;98%</td>
  </tr>""".format(color=color, num_favored=data["num_favored"], num_residues=data["num_residues"], favored_percentage=data["favored_percentage"]))

    if "rama_z" in self.summaries:
      #not yet available in json
      #work w/Oleg to get compatible output
      pass

    #MolProbity Score
    if "clashscore" in self.summaries and "ramalyze" in self.summaries and "rotalyze" in self.summaries:
      rows += 1
      from mmtbx.validation.utils import molprobity_score
      from mmtbx.validation.molprobity import percentile_mpscore
      mpscore = molprobity_score(self.summaries["clashscore"][model]["clashscore"],
                                 self.summaries["rotalyze"][model]["outlier_percentage"],
                                 self.summaries["ramalyze"][model]["favored_percentage"])
      pct_stats = percentile_mpscore.get_percentile_for_mpscore(mpscore, resolution=resolution)
      #pct_stats = {"minres":minres, "maxres":maxres, "count":nSamples, "percentile":pctRank}
      if pct_stats["percentile"] >= 66: color=green
      elif pct_stats["percentile"] < 33: color=red
      else: color=yellow
      range_text = "(N={count}, {minres}&Aring; - {maxres}&Aring;)".format(count=pct_stats["count"], minres=pct_stats["minres"], maxres=pct_stats["maxres"])
      protein_lines.append("""  <tr>
    <td>MolProbity score</td>
    <td colspan='2' bgcolor='{color}'>&nbsp;{mpscore:.2f}</td>
    <td>{percentile}<sup>th</sup> percentile {range_text}</td>
  </tr>""".format(color=color, mpscore=mpscore, percentile=pct_stats["percentile"], range_text=range_text))

    if "cbetadev" in self.summaries:
      rows += 1
      data = self.summaries["cbetadev"][model]
      if data["num_outliers"] == 0: color=green
      elif data["outlier_percentage"] >= 5: color=red
      else: color=yellow
      protein_lines.append("  <tr>")
      protein_lines.append("""    <td>C&beta; deviations &gt;0.25&Aring;</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers} / {num_cbeta_residues}</td>
    <td bgcolor='{color}'>&nbsp;{outlier_percentage:.2f}%</td>
    <td>Goal: 0</td>
  </tr>""".format(color=color, num_outliers=data["num_outliers"], num_cbeta_residues=data["num_cbeta_residues"], outlier_percentage=data["outlier_percentage"]))

    if "mp_bonds" in self.summaries and self.summaries["mp_bonds"][model]["num_total_protein"] > 0:
      rows += 1
      data = self.summaries["mp_bonds"][model]
      if data["num_total_protein"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_protein"]/data["num_total_protein"]*100.0
      if pct_outliers < 0.01: color=green
      elif pct_outliers >= 0.2: color=red
      else: color=yellow
      protein_lines.append("  <tr>")
      protein_lines.append("""    <td>Bad bonds</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers_protein} / {num_total_protein}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: 0%</td>
  </tr>""".format(color=color, num_outliers_protein=data["num_outliers_protein"], num_total_protein=data["num_total_protein"], pct_outliers=pct_outliers))

    if "mp_angles" in self.summaries and self.summaries["mp_angles"][model]["num_total_protein"] > 0:
      rows += 1
      data = self.summaries["mp_angles"][model]
      if data["num_total_protein"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_protein"]/data["num_total_protein"]*100.0
      if pct_outliers < 0.1: color=green
      elif pct_outliers >= 0.5: color=red
      else: color=yellow
      protein_lines.append("  <tr>")
      protein_lines.append("""    <td>Bad angles</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers_protein} / {num_total_protein}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: 0.1%</td>
  </tr>""".format(color=color, num_outliers_protein=data["num_outliers_protein"], num_total_protein=data["num_total_protein"], pct_outliers=pct_outliers))

    if rows == 0:
      return []
    #Overwrite the first <tr> with the Protein Geometry section cell
    protein_lines[0] = """  <tr>
    <td align='center' rowspan='{rows}'>Protein Geometry</td>""".format(rows=rows)
    return protein_lines

  def summary_table_peptide_omegas(self, model, show_all=False):
    data = self.summaries["omegalyze"][model]
    lines = []
    rows = 1
    if data["num_proline"] == 0:
      pct_cis_pro = 0.0
    else:
      pct_cis_pro = data["num_cis_proline"]/data["num_proline"]*100.0
    lines.append("""    <td>Cis Prolines</td>
    <td>&nbsp;{num_cis_proline} / {num_proline}</td>
    <td>&nbsp;{pct_cis_pro:.2f}</td>
    <td>Expected: &le;1 per chain, or &le;5%</td>
  </tr>""".format(num_cis_proline=data["num_cis_proline"], num_proline=data["num_proline"], pct_cis_pro=pct_cis_pro))

    if show_all or data["num_cis_general"] != 0:
      rows += 1
      if data["num_general"] == 0:
        pct_cis_general = 0.0
      else:
        pct_cis_general = data["num_cis_general"]/data["num_general"]*100.0
      if pct_cis_general <= 0.05: color=green
      elif pct_cis_general > 0.1: color=red
      else: color=yellow
      lines.append("""  <tr>
    <td>Cis nonProlines</td>
    <td bgcolor='{color}'>&nbsp;{num_cis_general} / {num_general}</td>
    <td bgcolor='{color}'>&nbsp;{pct_cis_general}%</td>
    <td>Goal: &lt;0.05%</td>
  </tr>""".format(color=color, num_cis_general=data["num_cis_general"], num_general=data["num_general"], pct_cis_general=pct_cis_general))

    num_twisted = data["num_twisted_proline"] + data["num_twisted_general"]
    if show_all or num_twisted != 0:
      rows += 1
      num_res = data["num_general"]+data["num_proline"]
      if num_res == 0:
        pct_twisted = 0.0
      else:
        pct_twisted = num_twisted/num_res*100.0
      if num_twisted == 0: color=green
      elif pct_twisted > 0.1: color=red
      else: color=yellow
      lines.append("""  <tr>
    <td>Twisted peptides</td>
    <td bgcolor='{color}'>&nbsp;{num_twisted} / {num_res}</td>
    <td bgcolor='{color}'>&nbsp;{pct_twisted}%</td>
    <td>Goal: 0</td>
  </tr>""".format(color=color, num_twisted=num_twisted, num_res=num_res, pct_twisted=pct_twisted))
    omega_lines = ["  <tr>","    <td align='center' rowspan='{rows}'>Peptide Omegas</td>".format(rows=len(lines))] + lines
    return omega_lines

  def summary_table_nucleic_acids(self, model):
    lines=[]
    rows = 0
    if "rna_puckers" in self.summaries:
      rows += 1
      data = self.summaries["rna_puckers"][model]
      if data["num_residues"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers"]/data["num_residues"]*100.0
      if data["num_outliers"] == 0: color=green
      elif pct_outliers > 5: color=red
      else: color=yellow
      lines.append("  <tr>")
      lines.append("""    <td>Probably wrong sugar puckers</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: 0</td>
  </tr>""".format(color=color, num_outliers=data["num_outliers"], num_residues=data["num_residues"], pct_outliers=pct_outliers))

    if "rna_suites" in self.summaries:
      rows += 1
      data = self.summaries["rna_suites"][model]
      if data["num_suites"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers"]/data["num_suites"]*100.0
      if pct_outliers <= 5: color=green
      elif pct_outliers > 15: color=red
      else: color=yellow
      lines.append("  <tr>")
      lines.append("""    <td>Bad backbone conformations</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers} / {num_suites}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: &le; 5%</td>
  </tr>""".format(color=color, num_outliers=data["num_outliers"], num_suites=data["num_suites"], pct_outliers=pct_outliers))

    #if "rna_bonds" in self.summaries:
    if "mp_bonds" in self.summaries and self.summaries["mp_bonds"][model]["num_total_na"] > 0:
      rows += 1
      data = self.summaries["mp_bonds"][model]
      if data["num_total_na"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_na"]/data["num_total_na"]*100.0
      if pct_outliers < 0.01: color=green
      elif pct_outliers >= 0.2: color=red
      else: color=yellow
      lines.append("  <tr>")
      lines.append("""    <td>Bad bonds</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers_na} / {num_total_na}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: 0%</td>
  </tr>""".format(color=color, num_outliers_na=data["num_outliers_na"], num_total_na=data["num_total_na"], pct_outliers=pct_outliers))

    if "mp_angles" in self.summaries and self.summaries["mp_angles"][model]["num_total_na"] > 0:
      rows += 1
      data = self.summaries["mp_angles"][model]
      if data["num_total_na"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_na"]/data["num_total_na"]*100.0
      if pct_outliers < 0.1: color=green
      elif pct_outliers >= 0.5: color=red
      else: color=yellow
      lines.append("  <tr>")
      lines.append("""    <td>Bad bonds</td>
    <td bgcolor='{color}'>&nbsp;{num_outliers_na} / {num_total_na}</td>
    <td bgcolor='{color}'>&nbsp;{pct_outliers:.2f}%</td>
    <td>Goal: &lt;0.1%</td>
  </tr>""".format(color=color, num_outliers_na=data["num_outliers_na"], num_total_na=data["num_total_na"], pct_outliers=pct_outliers))

    if rows == 0:
      return []
    lines[0] = """  <tr>
    <td align='center' rowspan={rows}>Nucleic Acid<br>Geometry</td>""".format(rows=rows)
    return lines

  def summary_table_low_res(self, model): #Right now, this just means 2 rows of CaBLAM
    if "cablam" not in self.summaries:
      return []
    lines = []
    lines.append("""  <tr>
    <td align='center' rowspan='2'>Low-resolution Criteria</td>""")
    data = self.summaries["cablam"][model]
    if data["cablam_outliers_percentage"] <=1: color=green
    elif data["cablam_outliers_percentage"] > 5: color=red
    else: color=yellow
    lines.append("""    <td>CaBLAM outliers</td>
    <td bgcolor='{color}'>&nbsp;{num_cablam_outliers} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{cablam_outliers_percentage:.1f}%</td>
    <td>Goal: &lt;1.0%</td></tr>""".format(
    color=color, num_cablam_outliers=data["num_cablam_outliers"], num_residues=data["num_residues"], cablam_outliers_percentage=data["cablam_outliers_percentage"]))

    if data["ca_geom_outliers_percentage"] <=0.5: color=green
    elif data["ca_geom_outliers_percentage"] > 1: color=red
    else: color=yellow
    lines.append("""  <tr>
    <td>CA Geometry outliers</td>
    <td bgcolor='{color}'>&nbsp;{num_ca_geom_outliers} / {num_residues}</td>
    <td bgcolor='{color}'>&nbsp;{ca_geom_outliers_percentage:.1f}%</td>
    <td>Goal: &lt;1.0%</td></tr>""".format(
    color=color, num_ca_geom_outliers=data["num_ca_geom_outliers"], num_residues=data["num_residues"], ca_geom_outliers_percentage=data["ca_geom_outliers_percentage"]))

    return lines

  def summary_table_additional(self, model):
    first_row=True
    lines = []
    rows=0
    if "chirals" in self.summaries:
      rows+=1
      #nothing to load yet, add later
      lines.append("  <tr>")
      lines.append()
    if "undowser" in self.summaries:
      rows+=1
      data = self.summaries["undowser"][model]
      if data["num_waters"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers"]/data["num_waters"]*100.0
      lines.append("  <tr>")
      lines.append("""    <td>Waters with clashes</td>
    <td>&nbsp;{num_outliers} / {num_waters}</td>
    <td>&nbsp;{pct_outliers:.2f}%</td>
    <td>See UnDowser table for details</td>
  </tr>""".format(
       num_outliers=data["num_outliers"], num_waters=data["num_waters"], pct_outliers=pct_outliers))
    if rows == 0:
      return []
    lines[0] = """  <tr>
    <td align='center' rowspan='{rows}'>Additional validations</td>""".format(rows=rows)
    return lines

  def clash_header(self, model, bootstrap=False):
    clashscore_val = self.summaries["clashscore"][model]["clashscore"]
    sortable_attr = " data-sortable='true'" if bootstrap else ""
    
    return (
        f"    <th{sortable_attr}>Clash &gt; 0.4&Aring;<br><small>Clashscore: "
        f"{clashscore_val:.2f}</small></th>"
    )

  def get_clashscore_header_content(self, model):
    clashscore_val = self.summaries["clashscore"][model]["clashscore"]
    return (
        f"Clash &gt; 0.4&Aring;<br><small>Clashscore: "
        f"{clashscore_val:.2f}</small>"
    )

  def rama_header(self, model, bootstrap=False):
    if bootstrap: 
      return "    <th data-sortable='true'>Ramachandran<br><small>Outliers: {num_outliers} out of {num_residues}</small></th>".format(
        num_outliers=self.summaries["ramalyze"][model]["num_outliers"], num_residues=self.summaries["ramalyze"][model]["num_residues"])
    else:
      return "    <th>Ramachandran<br><small>Outliers: {num_outliers} out of {num_residues}</small></th>".format(
        num_outliers=self.summaries["ramalyze"][model]["num_outliers"], num_residues=self.summaries["ramalyze"][model]["num_residues"])

  def get_ramalyze_header_content(self, model):
    summary = self.summaries["ramalyze"][model]
    return (
        f"Ramachandran<br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_residues']}</small>"
    )

  def rota_header(self, model, bootstrap=False):
    summary = self.summaries["rotalyze"][model]
    sortable_attr = " data-sortable='true'" if bootstrap else ""
    
    return (
        f"    <th{sortable_attr}>Rotamer<br><small>Poor rotamers: "
        f"{summary['num_outliers']} out of {summary['num_residues']}</small></th>"
    )

  def get_rotalyze_header_content(self, model):
    summary = self.summaries["rotalyze"][model]    
    return (
        f"Rotamer<br><small>Poor rotamers: "
        f"{summary['num_outliers']} out of {summary['num_residues']}</small>"
    )

  def cbdev_header(self, model, bootstrap=False):
    if bootstrap:
      return "   <th data-sortable='true'>C&beta; deviation<br><small>Outliers: {num_outliers} out of {num_cbeta_residues}</small></th>".format(
        num_outliers=self.summaries["cbetadev"][model]["num_outliers"], num_cbeta_residues=self.summaries["cbetadev"][model]["num_cbeta_residues"])
    else:
      return "   <th>C&beta; deviation<br><small>Outliers: {num_outliers} out of {num_cbeta_residues}</small></th>".format(
        num_outliers=self.summaries["cbetadev"][model]["num_outliers"], num_cbeta_residues=self.summaries["cbetadev"][model]["num_cbeta_residues"])
    
  def cablam_header(self, model, bootstrap=False):
    summary = self.summaries["cablam"][model]
    sortable_attr = " data-sortable='true'" if bootstrap else ""

    return (
        f"    <th{sortable_attr}>CaBLAM<br><small>Outliers: "
        f"{summary['num_cablam_outliers']} out of {summary['num_residues']}</small></th>"
    )

  def get_cablam_header_content(self, model):
    summary = self.summaries["cablam"][model]
    return (
        f"CaBLAM<br><small>Outliers: "
        f"{summary['num_cablam_outliers']} out of {summary['num_residues']}</small>"
    )
    
  def pperp_header(self, model, bootstrap=False):
    return "    <th>Base-P perp dist.<br><small>Outliers: {num_outliers} out of {num_residues}</small></th>".format(
      num_outliers=self.summaries["rna_puckers"][model]["num_outliers"], num_residues=self.summaries["rna_puckers"][model]["num_residues"])

  def suite_header(self, model, bootstrap=False):
    return "    <th>RNA suite conf.<br><small>Outliers: {num_outliers} out of {num_suites}</small></th>".format(
      num_outliers=self.summaries["rna_suites"][model]["num_outliers"], num_suites=self.summaries["rna_suites"][model]["num_suites"])

  def bonds_header(self, model, bootstrap=False):
    if bootstrap:
      return "    <th data-sortable='true'>Bond lengths<br><small>Outliers: {num_outliers} out of {num_total}</small></small></th>".format(
        num_outliers=self.summaries["mp_bonds"][model]["num_outliers"], num_total=self.summaries["mp_bonds"][model]["num_total"])
    else:
      return "    <th>Bond lengths<br><small>Outliers: {num_outliers} out of {num_total}</small></small></th>".format(
        num_outliers=self.summaries["mp_bonds"][model]["num_outliers"], num_total=self.summaries["mp_bonds"][model]["num_total"])

  def get_mp_bonds_header_content(self, model):
    summary = self.summaries["mp_bonds"][model]
    return (
        f"Bond lengths<br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_total']}</small>"
    )

  def angles_header(self, model, bootstrap=False):
    if bootstrap:
      return "    <th data-sortable='true'>Bond angles<br><small>Outliers: {num_outliers} out of {num_total}</small></th>".format(
        num_outliers=self.summaries["mp_angles"][model]["num_outliers"], num_total=self.summaries["mp_angles"][model]["num_total"])
    else:
      return "    <th>Bond angles<br><small>Outliers: {num_outliers} out of {num_total}</small></th>".format(
        num_outliers=self.summaries["mp_angles"][model]["num_outliers"], num_total=self.summaries["mp_angles"][model]["num_total"])
    
  def get_mp_angles_header_content(self, model):
    summary = self.summaries["mp_angles"][model]
    return (
        f"Bond angles<br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_total']}</small>"
    )

  def omegalyze_header(self, model, bootstrap=False):
    total_res = self.summaries["omegalyze"][model]["num_proline"] + self.summaries["omegalyze"][model]["num_general"]
    total_nontrans = self.summaries["omegalyze"][model]["num_cis_proline"] + self.summaries["omegalyze"][model]["num_twisted_proline"] + self.summaries["omegalyze"][model]["num_cis_general"] + self.summaries["omegalyze"][model]["num_twisted_general"]
    return "    <th>Cis Peptides<br><small>Non-Trans: {total_nontrans} out of {total_res}</small></th>".format(
      total_nontrans=total_nontrans, total_res=total_res)

  def make_multicrit_table(self, model, outliers_only=False):
    lines = []
    #lines.append(self.html_header())
    if model.strip():
      lines.append("<br><b>Multi-criterion table: Model {model}</b>".format(model))
    else:
      lines.append("<br><b>Multi-criterion table</b>")
    lines.append("<table width='100%' cellspacing='1' border='1'>")
    lines.append("  <tr align='center' bgcolor='#9999cc'>")
    #lines.append("    <th><small>Model</small></th>")

    lines.append("    <th>Chain</th>")
    lines.append("    <th>#</th>")
    lines.append("    <th>ins</th>")
    lines.append("    <th>Res</th>")
    lines.append("    <th>alt</th>")
    #for val_type in self.validation_types:
    table_order = self.get_table_order()
    for val_type in table_order:
      #lines.append("    <td>%s</td>" % val_type) #this needs more details
      #lines.append("    <th style='position:sticky; top:10px; background-color:#9999cc'>%s</th>" % val_type) #this needs more details
      if val_type == "clashscore":
        lines.append(self.clash_header(model))
      elif val_type == "ramalyze":
        lines.append(self.rama_header(model))
      elif val_type == "rotalyze":
        lines.append(self.rota_header(model))
      elif val_type == "cbetadev":
        lines.append(self.cbdev_header(model))
      elif val_type == "cablam":
        lines.append(self.cablam_header(model))
      elif val_type == "rna_puckers":
        lines.append(self.pperp_header(model))
      elif val_type == "rna_suites":
        lines.append(self.suite_header(model))
      elif val_type == "mp_bonds":
        lines.append(self.bonds_header(model))
      elif val_type == "mp_angles":
        lines.append(self.angles_header(model))
      elif val_type == "omegalyze":
        lines.append(self.omegalyze_header(model))
      else:
        lines.append("    <th>%s</th>" % val_type) #this needs more details
    lines.append("  </tr>")

    line_colors = ['#ffffff','#f0f0f0']
    count = 0

    for reskey in self.indices:
      if outliers_only and not self.data[reskey].has_outlier:
        continue
      color = line_colors[count % 2]
      lines.append(self.data[reskey].residue_cells(color))
      for alt in self.data[reskey].alternates:
        if alt != self.data[reskey].alternates[0]:
          lines.append("  <tr align='center' bgcolor='{color}'>".format(color=color))
        lines.append("    <td>{resname}</td>".format(resname=self.data[reskey].find_resname(alt=alt)))
        if alt in self.data[reskey].modeled_alternates:
          lines.append("    <td>{alt}</td>".format(alt=alt))
        else:
          lines.append("    <td>({alt})</td>".format(alt=alt))
        ##for val_type in self.validation_types:
        for val_type in table_order:
          if val_type == "clashscore":
            lines.append(self.data[reskey].clash_cell(alt, reskey))
          elif val_type == "ramalyze":
            lines.append(self.data[reskey].ramalyze_cell(alt))
          elif val_type == "rotalyze":
            lines.append(self.data[reskey].rotalyze_cell(alt))
          elif val_type == "cbetadev":
            lines.append(self.data[reskey].cbetadev_cell(alt))
          elif val_type == "cablam":
            lines.append(self.data[reskey].cablam_cell(alt))
          elif val_type == "rna_puckers":
            lines.append(self.data[reskey].base_p_perp_cell(alt))
          elif val_type == "rna_suites":
            lines.append(self.data[reskey].suitename_cell(alt))
          elif val_type == "mp_bonds":
            lines.append(self.data[reskey].bond_length_cell(alt))
          elif val_type == "mp_angles":
            lines.append(self.data[reskey].bond_angle_cell(alt))
          elif val_type == "omegalyze":
            lines.append(self.data[reskey].omegalyze_cell(alt))

        lines.append("  </tr>")
      count+=1
    lines.append("</table>")
    return "\n".join(lines)

  def make_multicrit_bootstrap_table(self, model, outliers_only=False):
    lines = []
    #lines.append(self.html_header())
    if model.strip():
      lines.append("<br><b>Multi-criterion table: Model {model}</b>".format(model))
    else:
      lines.append("<br><b>Multi-criterion table</b>")
    lines.append("<table id='multicrit-table' data-sort-empty-last='true' class='table-striped' width='100%' cellspacing='1' border='1' data-toggle='table'>")
    lines.append("  <thead>")
    lines.append("  <tr>")
    #lines.append("    <th><small>Model</small></th>")

    lines.append("    <th data-sortable='true'>Chain</th>")
    lines.append("    <th data-sortable='true'>#</th>")
    lines.append("    <th data-sortable='true'>ins</th>")
    lines.append("    <th data-sortable='true'>Res</th>")
    lines.append("    <th data-sortable='true'>alt</th>")
    #for val_type in self.validation_types:
    table_order = self.get_table_order()
    for val_type in table_order:
      #lines.append("    <td>%s</td>" % val_type) #this needs more details
      #lines.append("    <th style='position:sticky; top:10px; background-color:#9999cc'>%s</th>" % val_type) #this needs more details
      if val_type == "clashscore":
        lines.append(self.clash_header(model, bootstrap=True))
      elif val_type == "ramalyze":
        lines.append(self.rama_header(model, bootstrap=True))
      elif val_type == "rotalyze":
        lines.append(self.rota_header(model, bootstrap=True))
      elif val_type == "cbetadev":
        lines.append(self.cbdev_header(model, bootstrap=True))
      elif val_type == "cablam":
        lines.append(self.cablam_header(model, bootstrap=True))
      elif val_type == "rna_puckers":
        lines.append(self.pperp_header(model, bootstrap=True))
      elif val_type == "rna_suites":
        lines.append(self.suite_header(model, bootstrap=True))
      elif val_type == "mp_bonds":
        lines.append(self.bonds_header(model, bootstrap=True))
      elif val_type == "mp_angles":
        lines.append(self.angles_header(model, bootstrap=True))
      elif val_type == "omegalyze":
        lines.append(self.omegalyze_header(model, bootstrap=True))
      else:
        lines.append("    <th>%s</th>" % val_type) #this needs more details
    lines.append("  </tr>")
    lines.append("</thead>")

    lines.append("<tbody>")
    for reskey in self.indices:
      self.data[reskey].get_alternates() # make sure all residues have all alts needed for chart
    for reskey in self.indices:
      if outliers_only and not self.data[reskey].has_outlier:
        continue
      lines.append(self.data[reskey].as_html_row())
    lines.append("</tbody>")
    lines.append("</table>")
    return "\n".join(lines)


  def html_header(self):
    return """<!DOCTYPE html>
 <html lang="en">

  <head>
  <title>MolProbity Table Test</title>
  <style>
  th {position:sticky; top:0px; background-color:#9999cc}
  </style>
  </head>
  <body>
"""

  def html_header_d(self):
    return """<!DOCTYPE html>

<html lang="en">

  <head>

    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">

    <title>MolProbity Table Test</title>


    <!-- CSS
    –––––––––––––––––––––––––––––––––––––––––––––––––– -->
    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
    <!-- Custom CSS -->
    <link rel="stylesheet" href="/static/css/mystyle.css" type="text/css">

    <!-- JAVASCRIPT
    –––––––––––––––––––––––––––––––––––––––––––––––––– -->

    <!-- jQuery library -->
    <script
      src="https://code.jquery.com/jquery-3.4.1.min.js"
      integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo="
      crossorigin="anonymous">
    </script>
    <!-- Popper.JS -->
    <script
      src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js"
      integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49"
      crossorigin="anonymous">
    </script>

    <!-- Bootstrap JS -->
    <script
      src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"
      integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM"
      crossorigin="anonymous">
    </script>

    <!-- Maphighlight -->
    <script type="text/javascript" src="/static/js/jquery.maphilight.min.js"></script>
    <!-- Custom javascript -->
    <script
        type="text/javascript"
        src="/static/django_tables2_column_shifter/js/django_tables2_column_shifter.min.js">
    </script>
    <script type="text/javascript" src="/static/js/myjavascript.js"></script>

  </head>
  """

class residue():
  def __init__(self, rg, model_id, chain_id, resseq, icode):
    self.reskey = make_reskey(model_id, chain_id, resseq, icode)
    self.model_id = model_id
    self.chain_id = chain_id
    self.resseq = resseq
    self.icode = icode
    self.resnames = self.get_resnames_from_hierarchy(rg)
    self.resname = self.resnames[sorted(self.resnames.keys())[0]] #returns resname of first altloc alphabetically
    self.has_outlier = False
    self.validations = {}
    self.alternates = ['']
    self.modeled_alternates = self.get_alts_from_hierarchy(rg)
    self.colors = {"severe":"#ee4d4d",
                   "outlier":"#ff76a9",
                   "mild":"#ffb3cc"}

  def get_alternates(self):
    allalts = []+self.modeled_alternates
    for validation in self.validations:
      allalts += self.validations[validation]
    if not allalts:
      allalts.append('')
    self.alternates = sorted(set(allalts))
    #print(self.reskey, self.alternates)

  def get_resnames_from_hierarchy(self, rg):
    #Intended to handle the rare case where a residue has alternates with different residue names
    #May eventually need work to tie residue names to altlocs
    #This is assumed rare enough to be low priority
    resnames = {}
    for ag in rg.atom_groups():
      resnames[ag.altloc] = ag.resname
    return resnames

  def get_alts_from_hierarchy(self, rg):
    #Help track whether an alt came from the model or from calculations
    modeled_alternates = []
    for ag in rg.atom_groups():
      modeled_alternates.append(ag.altloc.strip())
    return modeled_alternates

  def find_resname(self, alt=''):
    res = self.resnames.get(alt)
    if res:
      return res
    else:
      #return self.resname+"*" #this version marks the guess/default
      return self.resname #default based on first alphabetical altloc

  def residue_cells(self,color=None):
    #residue identifiers, plus row management for alternates
    lines = []
    if color is None:
      lines.append("  <tr>")
    else:
      lines.append("  <tr align=center bgcolor='{color}'>".format(color=color))
    #lines.append("    <td rowspan={rows}>{resid}</td>".format(rows=len(self.alternates), resid=self.model_id))
    lines.append("    <td rowspan={rows}>{resid}</td>".format(rows=len(self.alternates), resid=self.chain_id))
    lines.append("    <td rowspan={rows}>{resid}</td>".format(rows=len(self.alternates), resid=self.resseq))
    lines.append("    <td rowspan={rows}>{resid}</td>".format(rows=len(self.alternates), resid=self.icode))

    return("\n".join(lines))

  def clash_cell(self, alt, reskey, bootstrap=False):
    clash_data_res = self.validations.get("clashscore")
    if clash_data_res is None:
      return "    <td>-</td>"
    clash_data = clash_data_res.get(alt)
    if clash_data is None:
      return "    <td>-</td>"
    clash_size = abs(clash_data[0]["overlap"])
    if clash_size < 0.5: color = self.colors["mild"] #implicitely clash_size>=0.4
    elif clash_size > 0.9: color = self.colors["severe"]
    else: color = self.colors["outlier"]
    line1 = "<td bgcolor={color}>{clash_size:.2f}&Aring;".format(color=color, clash_size=clash_size)
    if reskey == make_reskey(clash_data[0]["model_id"], clash_data[0]["chain_id"], clash_data[0]["resseq"], clash_data[0]["icode"]):
      line2 = "<br><small>{src_atom} with {trg_chain} {trg_resseq}{trg_icode} {trg_resname} {trg_atom} {trg_altloc}</small>".format(src_atom=clash_data[0]["name"],
                                                                                                                 trg_chain=clash_data[0]["target_chain_id"],
                                                                                                                 trg_resseq=clash_data[0]["target_resseq"],
                                                                                                                 trg_icode=clash_data[0]["target_icode"],
                                                                                                                 trg_resname=clash_data[0]["target_resname"],
                                                                                                                 trg_atom=clash_data[0]["target_name"],
                                                                                                                 trg_altloc=clash_data[0]["target_altloc"])
    else: #swap source and target in the cell
      line2 = "<br><small>{src_atom} with {trg_chain} {trg_resseq}{trg_icode} {trg_resname} {trg_atom} {trg_altloc}</small>".format(src_atom=clash_data[0]["target_name"],
                                                                                                                 trg_chain=clash_data[0]["chain_id"],
                                                                                                                 trg_resseq=clash_data[0]["resseq"],
                                                                                                                 trg_icode=clash_data[0]["icode"],
                                                                                                                 trg_resname=clash_data[0]["resname"],
                                                                                                                 trg_atom=clash_data[0]["name"],
                                                                                                                 trg_altloc=clash_data[0]["altloc"])
    return line1+line2+"</td>"

  def bond_length_cell(self, alt, bootstrap=False):
    geom_data_res = self.validations.get("mp_bonds")
    if geom_data_res is None:
      return "    <td>-</td>"
    geom_data = geom_data_res.get(alt)
    if geom_data is None:
      return "    <td>-</td>"
    geom_score = abs(geom_data[0]["score"])
    if geom_score >= 10:
      if bootstrap: color = " data-status='severe'"
      else: color = " bgcolor={color}".format(color=self.colors["severe"])
    else: 
      if bootstrap: color = " data-status='outlier'"
      else: color = " bgcolor={color}".format(color=self.colors["outlier"]) #implicitly score>=4
    line1 = "<td{color}>{outlier_count} OUTLIERS(S)".format(color=color, outlier_count=len(geom_data))
    line2 = "<br><small>worst is {atom1}--{atom2}: {geom_score:.1f} &sigma;</small>".format(
      atom1=geom_data[0]["atoms_name"][0].strip(), atom2=geom_data[0]["atoms_name"][1].strip(), geom_score=geom_score)
    return line1+line2+"</td>"

  def bond_angle_cell(self, alt, bootstrap=False):
    geom_data_res = self.validations.get("mp_angles")
    if geom_data_res is None:
      return "    <td>-</td>"
    geom_data = geom_data_res.get(alt)
    if geom_data is None:
      return "    <td>-</td>"
    geom_score = abs(geom_data[0]["score"])
    if geom_score >= 10:
      if bootstrap: color = " data-status='severe'"
      else: color = " bgcolor={color}".format(color=self.colors["severe"])
    else: 
      if bootstrap: color = " data-status='outlier'"
      else: color = " bgcolor={color}".format(color=self.colors["outlier"]) #implicitly score>=4
    line1 = "<td{color}>{outlier_count} OUTLIERS(S)".format(color=color, outlier_count=len(geom_data))
    line2 = "<br><small>worst is {atom1}-{atom2}-{atom3}: {geom_score:.1f} &sigma;</small>".format(
      atom1=geom_data[0]["atoms_name"][0].strip(), atom2=geom_data[0]["atoms_name"][1].strip(), atom3=geom_data[0]["atoms_name"][2].strip(), geom_score=geom_score)
    return line1+line2+"</td>"

  def ramalyze_cell(self, alt, bootstrap=False):
  #<td bgcolor='#ffb3cc'>Allowed (0.36%)<br><small>General / -156.1,27.4</small></td>
    rama_data_res = self.validations.get("ramalyze")
    if rama_data_res is None:
      return "    <td>-</td>"
    rama_data = rama_data_res.get(alt)
    if rama_data is None:
      return "    <td>-</td>"
    #linestart = "    <td "
    #print("\n".join(rama_data.keys()))
    if rama_data["rama_type"] == "Allowed":
      if bootstrap:
        linestart = "    <td data-status='allowed'>"
      else:
        linestart = "    <td bgcolor='#ffb3cc'>"
    elif rama_data["rama_type"].lower() == "outlier":
      if bootstrap:
        linestart = "    <td data-status='outlier'>"
      else:
        linestart = "    <td bgcolor='#ff76a9'>"
    else:
      linestart = "<td>"
    line1 = "{rama_type} ({score:.2f}%)".format(rama_type=rama_data["rama_type"], score=rama_data["score"])
    line2 = "<br><small>{res_type}/ {phi:.1f},{psi:.1f}</small>".format(res_type=rama_data["res_type_label"], phi=rama_data["phi"], psi=rama_data["psi"])
    return linestart+line1+line2+"</td>"

  def rotalyze_cell(self, alt, bootstrap=False):
    #<td bgcolor='#ffb3cc'>Allowed (0.4%) <i>p0</i><br><small>chi angles: 41.2,280.7</small></td>
    rota_data_res = self.validations.get("rotalyze")
    if rota_data_res is None:
      return "    <td>-</td>"
    rota_data = rota_data_res.get(alt)
    if rota_data is None:
      return "    <td>-</td>"

    if rota_data["evaluation"] == "Allowed":
      if bootstrap: linestart = "    <td data-status='allowed'>"
      else: linestart = "    <td bgcolor='#ffb3cc'>"
    elif rota_data["evaluation"].lower() == "outlier":
      if bootstrap: linestart = "    <td data-status='outlier'>"
      else: linestart = "    <td bgcolor='#ff76a9'>"
    else:
      linestart = "    <td>"

    if rota_data["outlier"]:
      line1 = "{rota_eval} ({score:.2f}%)".format(rota_eval=rota_data["evaluation"], score=rota_data["score"])
    else:
      line1 = "{rota_eval} ({score:.2f}%) <i>{rotamer}</i> ".format(rota_eval=rota_data["evaluation"], score=rota_data["score"], rotamer=rota_data["rotamer_name"])

    chi_angle_list = []
    for chi in rota_data["chi_angles"]:
      if chi is not None:
        chi_angle_list.append("%.1f" % chi)
    line2 = "<br><small>chi angles: {chi_angles}</small>".format(chi_angles=",".join(chi_angle_list))
    return linestart+line1+line2+"</td>"

  def cbetadev_cell(self, alt, bootstrap=False):
    #<td bgcolor='#ff76a9'>0.26&Aring;</td>
    cbdev_data_res = self.validations.get("cbetadev")
    if cbdev_data_res is None:
      return "    <td>-</td>"
    cbdev_data = cbdev_data_res.get(alt)
    if cbdev_data is None:
      return "    <td>-</td>"
    if cbdev_data["outlier"]:
      line = "    <td bgcolor='#ff76a9'>{dev_dist:.2f}&Aring;</td>".format(dev_dist=cbdev_data["deviation"])
    else:
      line = "    <td>{dev_dist:.2f}&Aring;</td>".format(dev_dist=cbdev_data["deviation"])
    return line

  def cablam_cell(self, alt, bootstrap=False):
    #<td bgcolor='#ffb3cc'> CaBLAM Disfavored   (1.671%)<br><small>                 </small></td>
    cablam_data_res = self.validations.get("cablam")
    if cablam_data_res is None:
      return "    <td>-</td>"
    cablam_data = cablam_data_res.get(alt)
    if cablam_data is None:
      return "    <td>-</td>"

    score = cablam_data["scores"]["cablam"]
    if cablam_data['outlier_type'] == "CA Geom Outlier":
      line1 = "    <td bgcolor={color}>CA Geom Outlier ({score:.3f}%)".format(color=self.colors["severe"], score=score)
    elif cablam_data['outlier_type'] == "CaBLAM Outlier":
      line1 = "    <td bgcolor={color}>CaBLAM Outlier ({score:.3f}%)".format(color=self.colors["outlier"], score=score)
    elif cablam_data['outlier_type'] == "CaBLAM Disfavored":
      line1 = "    <td bgcolor={color}>CaBLAM Disfavored ({score:.3f}%)".format(color=self.colors["mild"], score=score)
    else:
      line1 = "    <td>Favored ({score:.3f}%)".format(score=score)

    if cablam_data["feedback"] == "try alpha helix":
      line2 = "<br><small>alpha helix</small>"
    elif cablam_data["feedback"] == "try beta sheet":
      line2 = "<br><small>beta sheet</small>"
    elif cablam_data["feedback"] == "try three-ten":
      line2 = "<br><small>three-ten</small>"
    else:
      line2 = ""
    return line1+line2+"</td>"

  def omegalyze_cell(self, alt, bootstrap=False):
    omega_data_res = self.validations.get("omegalyze")
    if omega_data_res is None:
      return "    <td>-</td>"
    omega_data = omega_data_res.get(alt)
    if omega_data is None:
      return "    <td>-</td>"

    #Cis Pro gets text, but no color
    #Cis nonPro gets outlier pink
    #Twisted Anything gets severe red
    if omega_data["omega_type"] == "Twisted":
      linestart = "    <td bgcolor={color}>".format(color=self.colors["severe"])
    elif omega_data["omega_type"] == "Cis":
      if omega_data["resname"] == "PRO":
        linestart = "    <td>"
      else:
        linestart = "    <td bgcolor={color}>".format(color=self.colors["outlier"])
    else:
      return "    <td>-</td>"

    if omega_data["resname"] == "PRO":
      line1 = "{omega_type} PRO".format(omega_type = omega_data["omega_type"])
    else:
      line1 = "{omega_type} nonPRO".format(omega_type = omega_data["omega_type"])

    line2 = "<br><small>omega= {omega:.2f}</small>".format(omega=omega_data["omega"])
    return linestart+line1+line2+"</td>"

  def suitename_cell(self, alt, bootstrap=False):
    suite_data_res = self.validations.get("rna_suites")
    if suite_data_res is None:
      return "    <td>-</td>"
    suite_data = suite_data_res.get(alt)
    if suite_data is None:
      return "    <td>-</td>"

    if suite_data["outlier"]:
      line1 = "    <td bgcolor={color}>OUTLIER".format(color=self.colors["outlier"])
      line2 = ""
    else:
      line1 = "    <td>conformer: {suitename}".format(suitename=suite_data["cluster"])
      line2 = ""

    if suite_data["bin"] == "trig": #i.e. triaged
      line2 = "<br><small>(triaged {reason})</small>".format(reason = suite_data["reason"])
    elif suite_data["bin"] == "inc": #i.e. incomplete
      line2 = "<br><small>(incomplete)</small>"
    else: #suite_data["bin"] is something like "33 t" or "23 p", representing the two puckers and the m/p/t-staggered gamma dihedral
      ddg = suite_data["bin"][0] + "'," + suite_data["bin"][1] + "'," + suite_data["bin"][3]
      if suite_data["outlier"]:
        line2 = "<br><small>&delta;-1,&delta;,&gamma;={suite_bin}</small>".format(suite_bin=ddg)
      else:
        line2 = "<br><small>&delta;-1,&delta;,&gamma;={suite_bin} ; suiteness={suiteness:.2f}</small>".format(suite_bin=ddg, suiteness=suite_data["suiteness"])
    return line1+line2+"</td>"

  def base_p_perp_cell(self, alt, bootstrap=False):
    pperp_data_res = self.validations.get("rna_puckers")
    if pperp_data_res is None:
      return "    <td>-</td>"
    pperp_data = pperp_data_res.get(alt)
    if pperp_data is None:
      return "    <td>-</td>"

    outlier_types = []
    if pperp_data["is_delta_outlier"]:
      outlier_types.append("&delta;")
    if pperp_data["is_epsilon_outlier"]:
      outlier_types.append("&epsilon;")
    outlier_type = " & ".join(outlier_types)

    #diagnosis should not be calculated here! Put it in the validation itself!
    diagnosis="unknown"
    #if pperp_data["is_delta_outlier"]
    #if($perpdist < 2.9) //2.9A is dist cutoff for C2' vs C3' endo pucker
    #          $probpucker = "C2'-endo";
    #        else
    #          $probpucker = "C3'-endo";

    line1 = "    <td bgcolor={color}>suspect sugar pucker - {outlier_type} outlier".format(color=self.colors["outlier"], outlier_type=outlier_type)
    line2 = "<br><small>(P-perp distance implies {diagnosis})</small>".format(diagnosis=diagnosis)
    return line1+line2+"</td>"

class residue_bootstrap():
  def __init__(self, rg, model_id, chain_id, resseq, icode):
    self.reskey = make_reskey(model_id, chain_id, resseq, icode)
    self.chain_id = chain_id
    self.resseq = resseq
    self.icode = icode
    self.resnames = self._get_resnames_from_hierarchy(rg)
    self.resname = self.resnames[sorted(self.resnames.keys())[0]]
    self.validations = {}
    self.has_outlier = False
    self.alternates = ['']
    self.modeled_alternates = self.get_alts_from_hierarchy(rg)

  def get_alternates(self):
    """
    Combines alternates from the model and from validation results.
    """
    # Start with a set of alternates found in the model structure
    all_alts = set(self.modeled_alternates)
    
    # Add alternates found in any of the validation results
    for validation_dict in self.validations.values():
        all_alts.update(validation_dict.keys())

    # Ensure there's at least one entry ('') for non-alternate residues
    if not all_alts:
        all_alts.add('')
    
    # Sort the final, unique list of alternates
    self.alternates = sorted(list(all_alts))

  def get_alts_from_hierarchy(self, rg):
    #Help track whether an alt came from the model or from calculations
    modeled_alternates = []
    for ag in rg.atom_groups():
      modeled_alternates.append(ag.altloc.strip())
    return modeled_alternates

  def _get_resnames_from_hierarchy(self, rg):
    resnames = {}
    for ag in rg.atom_groups():
      resnames[ag.altloc] = ag.resname
    return resnames

  def find_resname(self, alt=''):
    return self.resnames.get(alt, self.resname)

  def as_dict(self, table_order):
    """Returns the residue's data as a dictionary, not an HTML string."""
    data = {
        'chain': self.chain_id,
        'resseq': self.resseq,
        'icode': self.icode,
    }

    # --- Identifier Cells ---
    resname_divs = [f"<div>{self.find_resname(alt)}</div>" for alt in self.alternates]
    alt_divs = []
    for alt in self.alternates:
        alt_display = alt if alt else '-'
        if alt not in self.modeled_alternates and alt:
            alt_display = f"({alt})"
        alt_divs.append(f"<div>{alt_display}</div>")
    
    data['resname_html'] = "".join(resname_divs)
    data['alt_html'] = "".join(alt_divs)

    for val_type in table_order:
        worst_sort_value = None
        divs = []
        
        method_name = f"_get_{val_type}_div"
        # Check if a helper method exists for this validation type
        if hasattr(self, method_name):
            helper_method = getattr(self, method_name)
            
            # Loop through alternates to build the content for this single cell
            for alt in self.alternates:
                # Call the correct helper for the current val_type
                div, sort_val = helper_method(alt)
                divs.append(div)
                if sort_val is not None and (worst_sort_value is None or sort_val > worst_sort_value):
                    worst_sort_value = sort_val
        
        # Store the raw sort value and the combined HTML for this validation type
        data[f'{val_type}_sort_value'] = worst_sort_value
        data[f'{val_type}_html'] = "".join(divs)
            
    return data

  # --- New Main Method to Generate a Single Table Row ---
  def as_html_row(self):
    # This dictionary will hold the combined text for each validation cell
    # "clashscore", "ramalyze", "rotalyze", "cbetadev", "cablam", "rna_puckers", 
    # "rna_suites", "mp_bonds", "mp_angles", "omegalyze"
    cell_contents = {key: [] for key in MASTER_ORDER}
    cell_sort_values = {key: None for key in MASTER_ORDER} # Dictionary to track sort values
    cell_contents.update({"resname": [], "alt": []})

    for alt in self.alternates:
      cell_contents["resname"].append(f"<div>{self.find_resname(alt)}</div>")
      alt_display = alt if alt else '-'
      # Check if the alt came from the model geometry and is not the blank placeholder
      if alt not in self.modeled_alternates and alt:
          alt_display = f"({alt})"
      cell_contents["alt"].append(f"<div>{alt_display}</div>")

      # Get content and status for each validation
      for val_type in MASTER_ORDER:
        method_name = f"_get_{val_type}_div"
        # Check if a helper method exists for this validation type
        if hasattr(self, method_name):
            helper_method = getattr(self, method_name)
            # Call the found method and append its result
            div, sort_val = helper_method(alt)
            cell_contents[val_type].append(div)
            if sort_val is not None and (cell_sort_values[val_type] is None or sort_val > cell_sort_values[val_type]):
              cell_sort_values[val_type] = sort_val
        else:
            # Fallback for any missing helper methods
            #cell_contents[val_type].append("<div>N/A</div>")
            print("")

    # Assemble the final HTML row
    row_parts = ["<tr>"]
    row_parts.append(f"<td>{self.chain_id}</td>")
    row_parts.append(f"<td>{self.resseq}</td>")
    row_parts.append(f"<td>{self.icode}</td>")
    row_parts.append(f"<td>{''.join(cell_contents['resname'])}</td>")
    row_parts.append(f"<td>{''.join(cell_contents['alt'])}</td>")
    
    # Add validation cells with their overall status
    for key in MASTER_ORDER:
      if len(cell_contents[key])>0:
        content = ''.join(cell_contents.get(key, []))
        sort_value = cell_sort_values.get(key)
        sort_attr = f'data-sort-value="{sort_value if sort_value is not None else ""}"'
        row_parts.append(f'<td {sort_attr}>{content}</td>')

    row_parts.append("</tr>")
    return "\n".join(row_parts)

  def _update_worst_status(self, statuses, key, new_status):
    # Helper to determine the most severe status for a cell
    severity = {"favored": 0, "allowed": 1, "mild": 1, "outlier": 2, "severe": 3}
    current_status = statuses.get(key, "favored")
    if severity.get(new_status, -1) > severity.get(current_status, -1):
        statuses[key] = new_status
  
  # --- New Helper Methods to Get Divs ---

  def _get_clashscore_div(self, alt):
    data_list = self.validations.get("clashscore", {}).get(alt)
    if not data_list:
        return ("<div>-</div>", None)

    # Data is a list of clashes, sorted worst to best; use the first one
    worst_clash = data_list[0]
    clash_size = abs(worst_clash["overlap"])
    
    # Determine status
    status = ""
    if clash_size < 0.5:
        status = "mild"
    elif clash_size > 0.9:
        status = "severe"
    else:
        status = "outlier"
    
    status_attr = f'data-status="{status}"' if status else ''
    
    # Format the main content string
    content_string = f"{clash_size:.2f} Å"
    
    # Format the details in the <small> tag
    if self.reskey == make_reskey(worst_clash["model_id"], worst_clash["chain_id"], worst_clash["resseq"], worst_clash["icode"]):
        details = "<small>{src_atom} with {trg_chain} {trg_resseq}{trg_icode} {trg_resname} {trg_atom} {trg_altloc}</small>".format(
            src_atom=worst_clash["name"],
            trg_chain=worst_clash["target_chain_id"],
            trg_resseq=worst_clash["target_resseq"],
            trg_icode=worst_clash["target_icode"],
            trg_resname=worst_clash["target_resname"],
            trg_atom=worst_clash["target_name"],
            trg_altloc=worst_clash["target_altloc"])
    else:
        details = "<small>{src_atom} with {trg_chain} {trg_resseq}{trg_icode} {trg_resname} {trg_atom} {trg_altloc}</small>".format(
            src_atom=worst_clash["target_name"],
            trg_chain=worst_clash["chain_id"],
            trg_resseq=worst_clash["resseq"],
            trg_icode=worst_clash["icode"],
            trg_resname=worst_clash["resname"],
            trg_atom=worst_clash["name"],
            trg_altloc=worst_clash["altloc"])
            
    content_string += f"<br>{details}"
    html_div = f"<div {status_attr}>{content_string}</div>"
    return (html_div, clash_size)

  def _get_ramalyze_div(self, alt):
    data = self.validations.get("ramalyze", {}).get(alt)
    if not data: return ("<div>-</div>", None)
    
    status_map = {"OUTLIER": "outlier", "ALLOWED": "allowed", "FAVORED": "favored"}
    status = status_map.get(data["rama_type"].upper(), "")
    status_attr = f'data-status="{status}"' if status else ''
    
    content = "{rama_type} ({score:.2f}%)<br><small>{res_type}/ {phi:.1f},{psi:.1f}</small>".format(**data)
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, data["score"])

  def _get_rotalyze_div(self, alt):
    data = self.validations.get("rotalyze", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    status_map = {"OUTLIER": "outlier", "ALLOWED": "allowed", "FAVORED": "favored"}
    status = status_map.get(data["evaluation"].upper(), "")
    status_attr = f'data-status="{status}"' if status else ''
    
    chi_angles = ",".join(["%.1f" % c for c in data.get("chi_angles", []) if c is not None])
    content = "{evaluation} ({score:.2f}%)<br><small>chi angles: {chi_angles}</small>".format(
        evaluation=data["evaluation"], score=data["score"], chi_angles=chi_angles
    )
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, data["score"])

  def _get_mp_bonds_div(self, alt):
    data = self.validations.get("mp_bonds", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    score = abs(data[0]["score"])
    status = "severe" if score >= 10 else "outlier"
    status_attr = f'data-status="{status}"'
    
    content = f"{len(data)} OUTLIER(S)<br><small>worst is {data[0]['atoms_name'][0].strip()}--{data[0]['atoms_name'][1].strip()}: {score:.1f} &sigma;</small>"
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, len(data))
    
  def _get_mp_angles_div(self, alt):
    data = self.validations.get("mp_angles", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    score = abs(data[0]["score"])
    status = "severe" if score >= 10 else "outlier"
    status_attr = f'data-status="{status}"'
    
    content = f"{len(data)} OUTLIER(S)<br><small>worst is {'-'.join(a.strip() for a in data[0]['atoms_name'])}: {score:.1f} &sigma;</small>"
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, len(data))

  def _get_cablam_div(self, alt):
    data = self.validations.get("cablam", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    outlier_type = data.get("outlier_type", "")
    status = ""
    if "CA Geom Outlier" in outlier_type:
        status = "severe"
    elif "CaBLAM Outlier" in outlier_type:
        status = "outlier"
    elif "Disfavored" in outlier_type:
        status = "mild"
    else:
      outlier_type = "Favored"
    
    status_attr = f'data-status="{status}"' if status else ''
    
    score = data.get("scores", {}).get("cablam", 0.0)
    content_string = f'{outlier_type} ({score:.3f}%)'
    
    feedback = data.get("feedback", "")
    if feedback:
        # Clean up the feedback text and add it
        feedback_text = feedback.replace("try ", "")
        content_string += f"<br><small>{feedback_text}</small>"
    html_div = f"<div {status_attr}>{content_string}</div>"
    return (html_div, score)

def loadModel(filename):  # from suitename/suites.py
  from iotbx.data_manager import DataManager
  dm = DataManager()  # Initialize the DataManager and call it dm
  dm.set_overwrite(True)  # tell the DataManager to overwrite files with the same name
  model = dm.get_model(filename)
  return model
#model.get_hierarchy()

def read_jsons_from_files(file_list):
  #accepts list of file paths
  #returns dictionary of jsons
  validations = {}
  for file_path in file_list:
    if not os.path.isfile(file_path):
      continue
    valfile = open(file_path)
    validation = json.load(valfile)
    if validation.get("mp_bonds"):
      validations["mp_bonds"] = validation["mp_bonds"]
      validations["mp_angles"] = validation["mp_angles"]
    else:
      validations[validation["validation_type"]] = validation
    valfile.close()
  return validations

if __name__ == '__main__':
  pdbfilepath = sys.argv[1]
  model = loadModel(pdbfilepath)
  h = model.get_hierarchy()
  #merged = make_storage_structure_from_hierarchy(h)
  merged = merged_validations(h)

  file_list = sys.argv[2:]
  validations = read_jsons_from_files(file_list)

  merged.add_validations(validations.values())
  merged.add_summaries(validations.values())

  print(merged.html_header())
  print(merged.make_summary_table(""))
  print(merged.make_multicrit_table("", outliers_only=True))

