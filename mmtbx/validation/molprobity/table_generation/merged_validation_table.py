from .residue_validation_data import ResidueValidationData
from .utils import make_reskey
from .html_builder import HtmlBuilder
from mmtbx.validation.molprobity import percentile_clashscore


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

VALIDATION_CATEGORIES = {
    "clashscore": "General",
    "ramalyze": "Protein",
    "rotalyze": "Protein",
    "cbetadev": "Protein",
    "cablam": "Protein",
    "omegalyze": "Protein",
    "mp_bonds": "General",
    "mp_angles": "General",
    "rna_puckers": "RNA",
    "rna_suites": "RNA",
}

# used for checkboxes that turn columns off and on in multicrit chart
VALIDATION_HEADER_NAMES = {
    "clashscore": "Clashes",
    "ramalyze": "Ramachandran",
    "rotalyze": "Rotamer",
    "cbetadev": "Cβ deviation",
    "cablam": "CaBLAM",
    "mp_bonds": "Bond lengths",
    "mp_angles": "Bond angles",
    "omegalyze": "Cis peptides",
    "rna_puckers": "RNA Puckers",
    "rna_suites": "RNA Suites"
}

class MergedValidationTable:
  """
  Main class to manage validation data and generate HTML tables.
  Uses residue_validation_data.py.
  Used by html_builder.py.
  """
  def __init__(self, hierarchy):
    self.data = {}
    self.indices = [] #reskeys in the same order as hierarchy
    self.validation_types = []
    self.summaries = {}
    self.model_ids = []
    self.html_builder = HtmlBuilder() # Uses the HTML generator
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
          self.data[reskey] = ResidueValidationData(rg, model_id, chain_id, resseq, icode)

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

  def get_validation_categories(self):
    return VALIDATION_CATEGORIES

  def get_validation_header_names(self):
    return VALIDATION_HEADER_NAMES

  def get_summary_table_data(self, model, resolution=2.0):
    """
    Generates a list of dictionaries containing all summary data,
    ready to be converted to JSON.
    """
    table_data = []

    # --- Clashscore ---
    if "clashscore" in self.summaries:
        data = self.summaries["clashscore"][model]
        pct_stats = percentile_clashscore.get_percentile_for_clashscore(data["clashscore"], resolution=resolution)
        status = self.get_clashscore_summary_stoplight(model)
        table_data.append({
            "category": "All-Atom Contacts",
            "metric": "Clashscore",
            "value": f"{data['clashscore']:.2f}",
            "percentile_html": f"{pct_stats['percentile']}<sup>th</sup> percentile",
            "status": status,
            "note": "100th percentile is the best among structures of comparable resolution; 0th is the worst."
        })

    # --- Protein Geometry ---
    if "rotalyze" in self.summaries:
        data = self.summaries["rotalyze"][model]
        status = self.get_rotalyze_summary_stoplight(model) # Assumes you have this helper
        table_data.append({
            "category": "Protein Geometry",
            "metric": "Poor rotamers",
            "value": f"{data['num_outliers']} / {data['num_residues']}",
            "percent_html": f"{data['outlier_percentage']:.1f}%",
            "goal": "&lt;0.3%",
            "status": status,
            "note": "Some validations count different numbers of residues. Gly residues have no sidechains for rotamers."
        })
        status = self.get_rotalyze_favored_summary_stoplight(model)
        table_data.append({
            "category": "Protein Geometry",
            "metric": "Favored rotamers",
            "value": f"{data['num_favored']} / {data['num_residues']}",
            "percent_html": f"{data['favored_percentage']:.1f}%",
            "goal": "&gt;98%",
            "status": status
        })

    if "ramalyze" in self.summaries:
        data = self.summaries["ramalyze"][model]
        status = self.get_ramalyze_summary_stoplight(model)
        table_data.append({
            "category": "Protein Geometry",
            "metric": "Ramachandran outliers",
            "value": f"{data['num_outliers']} / {data['num_residues']}",
            "percent_html": f"{data['outlier_percentage']:.1f}%",
            "goal": "&lt;0.05%",
            "status": status,
            "note": "Some validations count different numbers of residues. Ramachandran's &phi; or &psi; is incomplete at chain ends."
        })
        status = self.get_ramalyze_favored_summary_stoplight(model)
        table_data.append({
            "category": "Protein Geometry",
            "metric": "Ramachandran favored",
            "value": f"{data['num_favored']} / {data['num_residues']}",
            "percent_html": f"{data['favored_percentage']:.1f}%",
            "goal": "&gt;98%",
            "status": status
        })

    #MolProbity Score
    if "clashscore" in self.summaries and "ramalyze" in self.summaries and "rotalyze" in self.summaries:
      from mmtbx.validation.utils import molprobity_score
      from mmtbx.validation.molprobity import percentile_mpscore
      mpscore = molprobity_score(self.summaries["clashscore"][model]["clashscore"],
                                 self.summaries["rotalyze"][model]["outlier_percentage"],
                                 self.summaries["ramalyze"][model]["favored_percentage"])
      pct_stats = percentile_mpscore.get_percentile_for_mpscore(mpscore, resolution=resolution)
      #pct_stats = {"minres":minres, "maxres":maxres, "count":nSamples, "percentile":pctRank}
      status=None
      if pct_stats["percentile"] >= 66: status=green
      elif pct_stats["percentile"] < 33: status=red
      else: status=yellow
      table_data.append({
          "category": "Protein Geometry",
          "metric": "MolProbity score",
          "value": f"{mpscore:.2f}",
          "percent_html": f"{pct_stats['percentile']}<sup>th</sup> percentile",
          "status": status,
          "note": "MolProbity score combines the clashscore, rotamer, and Ramachandran evaluations into a single score, normalized to be on the same scale as X-ray resolution.  Lower is better."
        })

    if "cbetadev" in self.summaries:
      data = self.summaries["cbetadev"][model]
      status = self.get_cbetadev_summary_stoplight(model)
      table_data.append({
        "category": "Protein Geometry",
        "metric": "C&beta; deviations &gt;0.25&Aring;",
        "value": f"{data['num_outliers']} / {data['num_cbeta_residues']}",
        "percent_html": f"{data['outlier_percentage']:.1f}%",
        "goal": "0",
        "status": status,
        "note": "Some validations count different numbers of residues. Not all residues have a C&beta; atom, and won't be counted in the total."
      })

    # --- Bad Bonds ---
    if "mp_bonds" in self.summaries and self.summaries["mp_bonds"][model].get("num_total_protein", 0) > 0:
        data = self.summaries["mp_bonds"][model]
        num_outliers = data.get("num_outliers_protein", 0)
        num_total = data.get("num_total_protein", 0)

        pct_outliers = (num_outliers / num_total * 100.0) if num_total > 0 else 0.0

        status = self.get_mp_bonds_summary_stoplight(model)

        table_data.append({
            "category": "Protein Geometry",
            "metric": "Bad bonds",
            "value": f"{num_outliers} / {num_total}",
            "percent_html": f"{pct_outliers:.1f}%",
            "goal": "0%",
            "status": status
        })

    # --- Bad Angles ---
    if "mp_angles" in self.summaries and self.summaries["mp_angles"][model].get("num_total_protein", 0) > 0:
        data = self.summaries["mp_angles"][model]
        num_outliers = data.get("num_outliers_protein", 0)
        num_total = data.get("num_total_protein", 0)

        pct_outliers = (num_outliers / num_total * 100.0) if num_total > 0 else 0.0

        status = self.get_mp_angles_summary_stoplight(model)

        table_data.append({
            "category": "Protein Geometry",
            "metric": "Bad angles",
            "value": f"{num_outliers} / {num_total}",
            "percent_html": f"{pct_outliers:.1f}%",
            "goal": "&lt;0.1%",
            "status": status
        })

    # --- Peptide Omegas ---
    if "omegalyze" in self.summaries:
        data = self.summaries['omegalyze'][model]

        # --- Cis Prolines (Informational) ---
        num_proline = data.get("num_proline", 0)
        num_cis_proline = data.get("num_cis_proline", 0)
        pct_cis_pro = (num_cis_proline / num_proline * 100.0) if num_proline > 0 else 0.0
        table_data.append({
            "category": "Peptide Omegas",
            "metric": "Cis Prolines",
            "value": f"{num_cis_proline} / {num_proline}",
            "percent_html": f"{pct_cis_pro:.1f}%",
            "goal": "&le;1 per chain, or &le;5%",
            "status": None # No stoplight color for this informational row
        })

        # --- Cis non-Prolines ---
        num_general = data.get("num_general", 0)
        num_cis_general = data.get("num_cis_general", 0)
        pct_cis_general = (num_cis_general / num_general * 100.0) if num_general > 0 else 0.0

        status_cis_general = self.get_cis_nonpro_summary_stoplight(model)

        table_data.append({
            "category": "Peptide Omegas",
            "metric": "Cis non-Prolines",
            "value": f"{num_cis_general} / {num_general}",
            "percent_html": f"{pct_cis_general:.1f}%",
            "goal": "&lt;0.05%",
            "status": status_cis_general
        })

        # --- Twisted Peptides ---
        num_twisted = data.get("num_twisted_proline", 0) + data.get("num_twisted_general", 0)
        total_res = num_general + num_proline
        pct_twisted = (num_twisted / total_res * 100.0) if total_res > 0 else 0.0

        status_twisted = self.get_twisted_summary_stoplight(model)

        table_data.append({
            "category": "Peptide Omegas",
            "metric": "Twisted peptides",
            "value": f"{num_twisted} / {total_res}",
            "percent_html": f"{pct_twisted:.1f}%",
            "goal": "0",
            "status": status_twisted
        })

            # --- Low-resolution Criteria ---
        if "cablam" in self.summaries:
            data = self.summaries["cablam"][model]

            # --- CaBLAM Outliers ---
            outlier_pct = data.get("cablam_outliers_percentage", 0.0)
            status_cablam = yellow # Default to yellow
            if outlier_pct <= 1.0:
                status_cablam = green
            elif outlier_pct > 5.0:
                status_cablam = red

            table_data.append({
                "category": "Low-resolution Criteria",
                "metric": "CaBLAM outliers",
                "value": f'{data.get("num_cablam_outliers", 0)} / {data.get("num_residues", 0)}',
                "percent_html": f"{outlier_pct:.1f}%",
                "goal": "&lt;1.0%",
                "status": status_cablam
            })

            # --- CA Geometry Outliers ---
            ca_geom_pct = data.get("ca_geom_outliers_percentage", 0.0)
            status_ca_geom = yellow # Default to yellow
            if ca_geom_pct <= 0.5:
                status_ca_geom = green
            elif ca_geom_pct > 1.0:
                status_ca_geom = red

            table_data.append({
                "category": "Low-resolution Criteria",
                "metric": "CA Geometry outliers",
                "value": f'{data.get("num_ca_geom_outliers", 0)} / {data.get("num_residues", 0)}',
                "percent_html": f"{ca_geom_pct:.1f}%",
                "goal": "&lt;0.5%",
                "status": status_ca_geom
            })
    # ... You would continue this pattern for all other summary metrics ...
    # (Favored Rotamers, Favored Ramachandran, MolProbity Score, etc.)

    return table_data

  def get_clashscore_summary_stoplight(self, model, resolution=None):
    from mmtbx.validation.molprobity import percentile_clashscore
    color = None
    if "clashscore" in self.summaries:
      data = self.summaries["clashscore"][model]
      pct_stats = percentile_clashscore.get_percentile_for_clashscore(data["clashscore"], resolution=resolution)
      #pct_stats = {"minres":minres, "maxres":maxres, "count":nSamples, "percentile":pctRank}
      if pct_stats["percentile"] >= 66: color=green
      elif pct_stats["percentile"] < 33: color=red
      else: color=yellow
    return color

  def get_ramalyze_summary_stoplight(self, model):
    color = None
    if "ramalyze" in self.summaries:
      data = self.summaries["ramalyze"][model]
      if data["outlier_percentage"] <= 0.05: color=green
      elif data["outlier_percentage"] > 0.5 and data["num_outliers"] >= 2: color=red
      else: color=yellow
    return color

  def get_ramalyze_favored_summary_stoplight(self, model):
    color = None
    if "ramalyze" in self.summaries:
      data = self.summaries["ramalyze"][model]
      if data["favored_percentage"] >= 98: color=green
      elif data["favored_percentage"] < 95: color=red
      else: color=yellow
    return color

  def get_rotalyze_summary_stoplight(self, model):
    color = None
    if "rotalyze" in self.summaries:
      data = self.summaries["rotalyze"][model]
      if data["outlier_percentage"] <= 0.3: color=green
      elif data["outlier_percentage"] > 1.5: color=red
      else: color=yellow
    return color

  def get_rotalyze_favored_summary_stoplight(self, model):
    color = None
    if "rotalyze" in self.summaries:
      data = self.summaries["rotalyze"][model]
      if data["num_residues"] == 0:
        color = None #none would be more accurate; I think this correctly shows no color in analysis_result.html
      else:
        favored_pct = data["num_favored"]/data["num_residues"]*100.0
        if favored_pct >= 98: color=green
        elif favored_pct < 95: color=red
        else: color=yellow
    return color

  def get_cbetadev_summary_stoplight(self, model):
    color = None
    if "cbetadev" in self.summaries:
      data = self.summaries["cbetadev"][model]
      if data["num_outliers"] == 0: color=green
      elif data["outlier_percentage"] >= 5: color=red
      else: color=yellow
    return color

  def get_mp_bonds_summary_stoplight(self, model):
    color = None
    if "mp_angles" in self.summaries:
      data = self.summaries["mp_bonds"][model]
      if data["num_total_protein"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_protein"]/data["num_total_protein"]*100.0
      if pct_outliers < 0.01: color=green
      elif pct_outliers >= 0.2: color=red
      else: color=yellow
    return color

  def get_mp_angles_summary_stoplight(self, model):
    color = None
    if "mp_angles" in self.summaries:
      data = self.summaries["mp_angles"][model]
      if data["num_total_protein"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers_protein"]/data["num_total_protein"]*100.0
      if pct_outliers < 0.1: color=green
      elif pct_outliers >= 0.5: color=red
      else: color=yellow
    return color

  def get_cablam_summary_stoplight(self, model):
    color = None
    if "cablam" in self.summaries:
      data = self.summaries["cablam"][model]
      if data["cablam_outliers_percentage"] <=1: color=green
      elif data["cablam_outliers_percentage"] > 5: color=red
      else: color=yellow
    return color

  def get_rna_puckers_summary_stoplight(self, model):
    color = None
    if "rna_puckers" in self.summaries:
      data = self.summaries["rna_puckers"][model]
      if data["num_residues"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers"]/data["num_residues"]*100.0
      if data["num_outliers"] == 0: color=green
      elif pct_outliers > 5: color=red
      else: color=yellow
    return color

  def get_rna_suites_summary_stoplight(self, model):
    color = None
    if "rna_suites" in self.summaries:
      data = self.summaries["rna_suites"][model]
      if data["num_suites"] == 0:
        pct_outliers = 0.0
      else:
        pct_outliers = data["num_outliers"]/data["num_suites"]*100.0
      if pct_outliers <= 5: color=green
      elif pct_outliers > 15: color=red
      else: color=yellow
    return color

  def get_cis_nonpro_summary_stoplight(self, model):
    if "omegalyze" in self.summaries:
      data = self.summaries["omegalyze"][model]
      num_general = data.get("num_general", 0)
      num_cis_general = data.get("num_cis_general", 0)
      pct_cis_general = (num_cis_general / num_general * 100.0) if num_general > 0 else 0.0

      color = yellow # Default to yellow
      if pct_cis_general <= 0.05:
        color = green
      elif pct_cis_general > 0.1:
        color = red
    return color

  def get_twisted_summary_stoplight(self, model):
    if "omegalyze" in self.summaries:
      data = self.summaries["omegalyze"][model]
      num_general = data.get("num_general", 0)
      num_proline = data.get("num_proline", 0)
      num_twisted = data.get("num_twisted_proline", 0) + data.get("num_twisted_general", 0)
      total_res = num_general + num_proline
      pct_twisted = (num_twisted / total_res * 100.0) if total_res > 0 else 0.0

      color = yellow # Default to yellow
      if num_twisted == 0:
          color = green
      elif pct_twisted > 0.1:
          color = red
    return color

  # returns the worst color out of all the omegalyze categories.
  # This is used for the color of the header of multicrit table.
  def get_omegalyze_summary_stoplight(self, model):
    if "omegalyze" not in self.summaries:
      return None

    data = self.summaries["omegalyze"][model]

    # Define severity levels for the statuses
    severity = {"green": 1, "yellow": 2, "red": 3}
    worst_severity = severity["green"] # Start with the best possible status

    # --- Check Cis non-Prolines ---
    num_general = data.get("num_general", 0)
    if num_general > 0:
      pct_cis_general = data.get("num_cis_general", 0) / num_general * 100.0

      current_severity = severity["green"]
      if pct_cis_general > 0.1:
        current_severity = severity["red"]
      elif pct_cis_general > 0.05:
        current_severity = severity["yellow"]

      if current_severity > worst_severity:
        worst_severity = current_severity
    # --- Check Twisted Peptides ---
    num_twisted = data.get("num_twisted_proline", 0) + data.get("num_twisted_general", 0)
    if num_twisted > 0:
      num_res = num_general + data.get("num_proline", 0)
      pct_twisted = num_twisted / num_res * 100.0 if num_res > 0 else 0.0

      current_severity = severity["green"]
      if pct_twisted > 0.1:
        current_severity = severity["red"]
      else:
        current_severity = severity["yellow"]

      if current_severity > worst_severity:
        worst_severity = current_severity
    # --- Return the color corresponding to the worst status found ---
    if worst_severity == severity["red"]:
      return red
    elif worst_severity == severity["yellow"]:
      return yellow
    else:
      return green

  def finalize_data(self):
    """Finalize alternates after all data is loaded."""
    for res in self.residues.values():
        res.get_alternates()

  def _build_residue_dict(self, residue_data, table_order):
    """
    Builds the final data dictionary for a single residue row,
    combining raw data (for sorting) and formatted HTML (for display).
    """
    # 1. Get the raw, unformatted data for the residue
    raw_data = residue_data.get_row_data(table_order)

    # 2. Use the HtmlBuilder to generate the formatted HTML parts
    # This returns a dictionary of {'html': "...", 'sort_value': 12.34}
    formatted_cells = self.html_builder.build_all_cells_for_residue(raw_data, table_order)

    # 3. Assemble the final dictionary for the JSON response
    final_dict = {
      'chain': raw_data['chain_id'],
      'resseq': raw_data['resseq'],
      'icode': raw_data['icode'],
      'resname_html': formatted_cells['resname']['html'],
      'alt_html': formatted_cells['alt']['html']
    }

    # Add the sort and html values for each validation type
    for val_type in table_order:
      final_dict[f'{val_type}_sort_value'] = formatted_cells[val_type]['sort_value']
      final_dict[f'{val_type}_html'] = formatted_cells[val_type]['html']

    return final_dict

  def get_all_header_data(self, model):
    """
    Dynamically orchestrates the creation of all header content.
    """
    header_data = {}
    table_order = self.get_table_order()

    for val_type in table_order:
      # 1. Construct the names of the required methods dynamically
      stoplight_method_name = f"get_{val_type}_summary_stoplight"
      header_content_method_name = f"get_{val_type}_header_content"

      # 2. Check that the required methods exist
      if not hasattr(self, stoplight_method_name):
        raise AttributeError(
          f"Missing required stoplight method: '{stoplight_method_name}' in MergedValidationTable"
        )
      if not hasattr(self.html_builder, header_content_method_name):
        raise AttributeError(
          f"Missing required header content method: '{header_content_method_name}' in HtmlBuilder"
        )

      # 3. Get the methods using getattr()
      stoplight_method = getattr(self, stoplight_method_name)
      header_content_method = getattr(self.html_builder, header_content_method_name)

      # 4. Execute the logic
      color = stoplight_method(model)
      header_html = header_content_method(self.summaries, model, color)

      header_data[val_type] = header_html

    return header_data

  def get_summary_html(self, model):
    """
    Generates the HTML for the summary table.
    """
    summary_data = self.get_summary_table_data(model)
    return self.html_builder.build_summary_table(summary_data)

  def get_multicrit_html(self, model):
    """
    Generates the HTML for the main multi-criterion table.
    """
    # 1. Get the dynamic header content
    header_data = self.get_all_header_data(model)

    # 2. Build the table header HTML
    header_html = self.html_builder.build_multicrit_header(header_data)

    # 3. Build the table body HTML by processing each residue
    table_order = self.get_table_order()
    body_rows = []
    for residue_data in self.data.values():
      body_rows.append(
        self.html_builder.build_residue_row(residue_data, table_order)
      )

    # 4. Combine header and body into the final table
    return (
      f'{header_html}'
      f'<tbody>{"".join(body_rows)}</tbody>'
    )
