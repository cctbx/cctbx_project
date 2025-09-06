from .utils import make_reskey

class HtmlBuilder:
  """Generates all HTML for the validation tables."""

  def _get_clashscore_div(self, validations, alt):
    data_list = validations.get("clashscore", {}).get(alt)
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
    reskey = validations['reskey']
    if reskey == make_reskey(worst_clash["model_id"], worst_clash["chain_id"], worst_clash["resseq"], worst_clash["icode"]):
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

  def _get_ramalyze_div(self, validations, alt):
    data = validations.get("ramalyze", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    status_map = {"OUTLIER": "outlier", "ALLOWED": "allowed", "FAVORED": "favored"}
    status = status_map.get(data["rama_type"].upper(), "")
    status_attr = f'data-status="{status}"' if status else ''

    content = "<span>{rama_type}<span class='cell-detail'> ({score:.2f}%)</span></span><span><small>{res_type}<span class='cell-detail'>: {phi:.1f}&phi;,{psi:.1f}&psi;</small></span></span>".format(**data)
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, data["score"])

  def _get_rotalyze_div(self, validations, alt):
    data = validations.get("rotalyze", {}).get(alt)
    if not data:
      return ("<div>-</div>", None)

    status_map = {"OUTLIER": "outlier", "ALLOWED": "allowed", "FAVORED": "favored"}
    status = status_map.get(data["evaluation"].upper(), "")
    status_attr = f'data-status="{status}"' if status else ''

    # Conditionally add the rotamer name in italics if it's not an outlier
    rotamer_display = ""
    if data["evaluation"].upper() != "OUTLIER":
      rotamer_name = data.get("rotamer_name", "")
      if rotamer_name:
        rotamer_display = f": <i>{rotamer_name}</i>"

    chi_angles = ", ".join([f"{round(c)}&deg;" for c in data.get("chi_angles", []) if c is not None])

    # Use an f-string for better readability
    content = (
      f"<span>{data['evaluation']}{rotamer_display}"
      f"<span class='cell-detail'> ({data['score']:.2f}%)</span></span>"
      f"<span class='cell-detail'><small>chi angles: {chi_angles}</small></span>"
    )

    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, data.get("score"))

  def _get_mp_bonds_div(self, validations, alt):
    cutoff = 4
    data = validations.get("mp_bonds", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    score = abs(data[0]["score"])
    if score < cutoff: return ("<div>-</div>", None)
    status = "severe" if score >= 10 else "outlier"
    status_attr = f'data-status="{status}"'
    outliers = sum(1 for item in data if item['outlier'])

    content = f"{outliers} OUTLIER(S)<br><small>worst is {data[0]['atoms_name'][0].strip()}--{data[0]['atoms_name'][1].strip()}: {score:.1f} &sigma;</small>"
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, len(data))

  def _get_mp_angles_div(self, validations, alt):
    cutoff = 4
    data = validations.get("mp_angles", {}).get(alt)
    if not data: return ("<div>-</div>", None)

    score = abs(data[0]["score"])
    if score < cutoff: return ("<div>-</div>", None)
    status = "severe" if score >= 10 else "outlier"
    status_attr = f'data-status="{status}"'
    outliers = sum(1 for item in data if item['outlier'])

    content = f"{outliers} OUTLIER(S)<br><small>worst is {'-'.join(a.strip() for a in data[0]['atoms_name'])}: {score:.1f} &sigma;</small>"
    html_div = f"<div {status_attr}>{content}</div>"
    return (html_div, len(data))

  def _get_cbetadev_div(self, validations, alt):
    """Generates the div for C-beta deviation analysis."""
    data = validations.get("cbetadev", {}).get(alt)
    if not data:
      return ("<div>-</div>", None)

    deviation = data.get("deviation")
    status = "outlier" if data.get("outlier") else ""
    status_attr = f'data-status="{status}"' if status else ''

    # Format the content, adding the Angstrom symbol
    content = f"{deviation:.2f}&Aring;"

    html_div = f"<div {status_attr}>{content}</div>"

    # Return the HTML and the raw deviation for sorting
    return (html_div, deviation)

  def _get_cablam_div(self, validations, alt):
    data = validations.get("cablam", {}).get(alt)
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
    content_string = f'<span>{outlier_type}<span class="cell-detail"> ({score:.3f}%)</span></span>'

    feedback = data.get("feedback", "")
    if feedback:
        # Clean up the feedback text and add it
        feedback_text = feedback.replace("try ", "")
        content_string += f"<small>{feedback_text}</small>"
    html_div = f"<div {status_attr}>{content_string}</div>"
    return (html_div, score)

  def _get_omegalyze_div(self, validations, alt):
    """Generates the div for Omegalyze analysis."""
    data = validations.get("omegalyze", {}).get(alt)
    if not data:
        return ("<div>-</div>", None)

    omega_type = data.get("omega_type")
    resname = data.get("resname")
    omega_value = data.get("omega", 0.0)

    status = ""
    sort_value = 0 # Lower is better, so Trans peptides will sort first.

    if omega_type == "Twisted":
        status = "severe"
        sort_value = 3 # Highest severity
    elif omega_type == "Cis":
        if resname != "PRO":
            status = "outlier"
            sort_value = 2 # Medium severity
        else: # Cis-proline is notable but not an outlier
            sort_value = 1 # Low severity

    status_attr = f'data-status="{status}"' if status else ''

    proline_label = "PRO" if resname == "PRO" else "nonPRO"
    content = (
        f"<span>{omega_type} {proline_label}"
        f"<span class='cell-detail'><br><small>omega={omega_value:.2f}&deg;</small></span></span>"
    )

    html_div = f"<div {status_attr}>{content}</div>"

    return (html_div, sort_value)

  def build_all_cells_for_residue(self, raw_residue_data, table_order):
    """
    Takes raw data for one residue and returns a dictionary containing
    the formatted HTML and raw sort values for each of its cells.
    """
    formatted_cells = {}

    # --- Identifier Cells ---
    resname_divs = []
    alt_divs = []
    for alt in raw_residue_data['alternates']:
      alt_display = alt if alt else '-'
      if alt not in raw_residue_data['modeled_alternates'] and alt:
          alt_display = f"({alt})"
      alt_divs.append(f"<div>{alt_display}</div>")

    formatted_cells['resname'] = {'html': "".join(resname_divs), 'sort_value': None}
    formatted_cells['alt'] = {'html': "".join(alt_divs), 'sort_value': None}

    # --- Validation Cells ---
    for val_type in table_order:
      worst_sort_value = None
      divs = []

      method_name = f"_get_{val_type}_div"
      # Check if a helper method exists for this validation type
      if hasattr(self, method_name):
        helper_method = getattr(self, method_name)

        # Loop through alternates to build the content for this single cell
        for alt in raw_residue_data['alternates']:
          # Call the correct helper for the current val_type
          div, sort_val = helper_method(raw_residue_data, alt)
          divs.append(div)
          if sort_val is not None and (worst_sort_value is None or sort_val > worst_sort_value):
            worst_sort_value = sort_val

      formatted_cells[val_type] = {
          'html': "".join(divs),
          'sort_value': worst_sort_value
      }

    return formatted_cells

  def build_residue_row(self, residue_row_data, table_order):
    """
    Builds a single <tr> HTML string from a dictionary of raw residue data.
    """
    # "clashscore", "ramalyze", "rotalyze", "cbetadev", "cablam", "rna_puckers",
    # "rna_suites", "mp_bonds", "mp_angles", "omegalyze"
    cell_contents = {key: [] for key in MASTER_ORDER}
    cell_sort_values = {key: None for key in MASTER_ORDER} # Dictionary to track sort values
    cell_contents.update({"resname": [], "alt": []})

    for alt in residue_row_data['alternates']:
      cell_contents["resname"].append(f"<div>{residue_row_data['resnames'].get(alt, 'UNK')}</div>")
      alt_display = alt if alt else '-'
      # Check if the alt came from the model geometry and is not the blank placeholder
      if alt not in residue_row_data['modeled_alternates'] and alt:
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

  def get_clashscore_header_content(self, summaries, model, stoplight_color):
    clashscore_val = summaries["clashscore"][model]["clashscore"]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Clash &gt; 0.4&Aring;<span class='cell-detail'><br><small>Clashscore: "
        f"{clashscore_val:.2f}</small></span></div>"
    )

  def get_ramalyze_header_content(self, summaries, model, stoplight_color):
    summary = summaries["ramalyze"][model]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Ramachandran<span class='cell-detail'><br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_residues']}</small></span></div>"
    )

  def get_rotalyze_header_content(self, summaries, model, stoplight_color):
    summary = summaries["rotalyze"][model]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Rotamer<span class='cell-detail'><br><small>Poor rotamers: "
        f"{summary['num_outliers']} out of {summary['num_residues']}</small></span></div>"
    )

  def get_cbetadev_header_content(self, summaries, model, stoplight_color):
    summary = summaries["cbetadev"][model]
    num_outliers = summary["num_outliers"]
    num_residues = summary["num_cbeta_residues"]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>C&beta; deviation<span class='cell-detail'><br><small>Outliers: "
        f"{num_outliers} out of {num_residues}</small></span></div>"
    )

  def get_cablam_header_content(self, summaries, model, stoplight_color):
    summary = summaries["cablam"][model]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>CaBLAM<span class='cell-detail'><br><small>Outliers: "
        f"{summary['num_cablam_outliers']} out of {summary['num_residues']}</small></span></div>"
    )

  def pperp_header(self, summaries, model, stoplight_color):
    return "    <th>Base-P perp dist.<br><small>Outliers: {num_outliers} out of {num_residues}</small></th>".format(
      num_outliers=summaries["rna_puckers"][model]["num_outliers"], num_residues=self.summaries["rna_puckers"][model]["num_residues"])

  def suite_header(self, summaries, model, stoplight_color):
    return "    <th>RNA suite conf.<br><small>Outliers: {num_outliers} out of {num_suites}</small></th>".format(
      num_outliers=self.summaries["rna_suites"][model]["num_outliers"], num_suites=self.summaries["rna_suites"][model]["num_suites"])

  def get_mp_bonds_header_content(self, summaries, model, stoplight_color):
    summary = summaries["mp_bonds"][model]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Bond lengths<span class='cell-detail'><br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_total']}</small></span></div>"
    )

  def get_mp_angles_header_content(self, summaries, model, stoplight_color):
    summary = summaries["mp_angles"][model]
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Bond angles<span class='cell-detail'><br><small>Outliers: "
        f"{summary['num_outliers']} out of {summary['num_total']}</small></span></div>"
    )

  def get_omegalyze_header_content(self, summaries, model, stoplight_color):
    summary = summaries["omegalyze"][model]

    # Calculate totals safely using .get() with a default of 0
    num_proline = summary.get("num_proline", 0)
    num_general = summary.get("num_general", 0)
    total_res = num_proline + num_general

    total_nontrans = (
        summary.get("num_cis_proline", 0) +
        summary.get("num_twisted_proline", 0) +
        summary.get("num_cis_general", 0) +
        summary.get("num_twisted_general", 0)
    )
    # Create the inner HTML content
    return (
        f"<div class='th-inner' data-status='{stoplight_color}'>Peptide Omegas<span class='cell-detail'><br><small>Non-Trans: "
        f"{total_nontrans} out of {total_res}</small></span></div>"
    )

  def build_all_headers(self, validation_table, table_order, model):
    """Builds a dictionary of all header content."""
    header_content = {}
    for val_type in table_order:
      method_name = f"get_{val_type}_header_content"
      if hasattr(self, method_name):
        # Call the method on this HtmlBuilder instance
        header_method = getattr(self, method_name)
        header_content[val_type] = header_method(validation_table.summaries, model)
    return header_content

  def build_summary_table(self, summary_data):
    """
    Builds the summary table HTML from a list of data dictionaries,
    with support for rowspan and popover help icons.
    """
    if not summary_data:
      return ""

    # 1. Group the flat data list by category
    grouped_data = {}
    for row in summary_data:
      category = row.get("category", "Other")
      if category not in grouped_data:
        grouped_data[category] = []
      grouped_data[category].append(row)

    # 2. Build the HTML from the grouped data
    summary_html_parts = [
      '<div class="result-box" style="margin-bottom: 1em;">',
      '<h4>Summary Statistics</h4>',
      '<table id="summary-table" class="table table-sm table-bordered">',
      '<thead><tr><th>Category</th><th>Metric</th><th>Value</th><th>Percent</th><th>Goal</th></tr></thead>',
      '<tbody>'
    ]

    # 3. Loop through each category to build the rows
    for category, rows in grouped_data.items():
      rowspan = len(rows)

      # --- Handle the FIRST row of the category ---
      first_row = rows[0]
      help_icon_html = self._create_help_icon_html(first_row)

      summary_html_parts.append('<tr>')
      summary_html_parts.append(f'<td rowspan="{rowspan}" style="vertical-align: middle;">{category}</td>')
      summary_html_parts.append(f"<td><div>{first_row.get('metric', '')}{help_icon_html}</div></td>")
      summary_html_parts.append(f'<td><div data-status="{first_row.get("status", "")}">{first_row.get("value", "")}</div></td>')
      summary_html_parts.append(f'<td><div data-status="{first_row.get("status", "")}">{first_row.get("percentile_html") or first_row.get("percent_html", "")}</div></td>')
      summary_html_parts.append(f"<td>{first_row.get('goal', '')}</td>")
      summary_html_parts.append('</tr>')

      # --- Handle SUBSEQUENT rows for this category ---
      for i in range(1, rowspan):
        subsequent_row = rows[i]
        help_icon_html = self._create_help_icon_html(subsequent_row)

        summary_html_parts.append('<tr>')
        summary_html_parts.append(f"<td><div>{subsequent_row.get('metric', '')}{help_icon_html}</div></td>")
        summary_html_parts.append(f'<td><div data-status="{subsequent_row.get("status", "")}">{subsequent_row.get("value", "")}</div></td>')
        summary_html_parts.append(f'<td><div data-status="{subsequent_row.get("status", "")}">{subsequent_row.get("percentile_html") or subsequent_row.get("percent_html", "")}</div></td>')
        summary_html_parts.append(f"<td>{subsequent_row.get('goal', '')}</td>")
        summary_html_parts.append('</tr>')

    summary_html_parts.extend(['</tbody>', '</table>', '</div>'])
    return "\n".join(summary_html_parts)

  def _create_help_icon_html(self, row_data):
    """Creates the HTML for a popover help icon if a note exists."""
    note = row_data.get("note")
    if note:
      return (
        f'<a tabindex="0" class="ml-2" role="button" data-toggle="popover" '
        f'data-trigger="focus" title="Note" data-content="{note}">(?)</a>'
      )
    return ""

  def build_multicrit_header(self, header_data):
    """
    Builds the <thead> section for the multi-criterion table.
    """
    header_cells = [
      "<th data-sortable='true'>Chain</th>",
      "<th data-sortable='true'>#</th>",
      "<th data-sortable='true'>ins</th>",
      "<th data-sortable='true'>Res</th>",
      "<th data-sortable='true'>alt</th>"
    ]

    # Add the dynamically generated header content
    for val_type in table_order:
      header_html_div = header_data.get(val_type, {}).get('html', f"<div>{val_type.title()}</div>")
      # You can add the data-sorter attribute here if needed
      header_cells.append(f"<th>{header_html_div}</th>")

    return f"<thead><tr>{''.join(header_cells)}</tr></thead>"
