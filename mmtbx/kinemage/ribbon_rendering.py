"""Ribbon rendering for kinemage output.

Extracted from mmtbx/programs/ribbons.py into standalone classes and functions
that can be reused by both the standalone ribbon program and the MolProbity
kinemage validation pipeline.

Uses scitbx.matrix.col throughout instead of numpy for consistency with the
rest of the kinemage code.
"""

from __future__ import absolute_import, division, print_function
from scitbx.matrix import col
from mmtbx.kinemage.nrubs import Triple, NRUBS
from mmtbx.kinemage.ribbons import (
    find_contiguous_protein_residues, find_contiguous_nucleic_acid_residues,
    make_protein_guidepoints, make_nucleic_acid_guidepoints,
    untwist_ribbon, swap_edge_and_face, _FindNamedAtomInResidue,
    _IsNucleicAcidResidue, chain_has_DNA, chain_has_RNA)


class Crayon(object):
  """Controls per-point coloring for ribbon kinemage output.

  Can produce rainbow coloring (smoothly varying across a chain) or a
  constant color string.
  """

  _rainbowColors = [
      "blue", "sky", "cyan", "sea", "green",
      "lime", "yellow", "gold", "orange", "red"]

  def __init__(self, do_rainbow, color=None):
    self._do_rainbow = do_rainbow
    self._color = color

  def for_ribbon(self, start_guide):
    """Update color based on the ribbon position (rainbow mode only)."""
    if self._do_rainbow:
      res = start_guide.prevRes
      chain = res.parent()
      first_resid = chain.residue_groups()[0].resseq_as_int()
      last_resid = chain.residue_groups()[-1].resseq_as_int()
      denom = last_resid - first_resid
      if denom == 0:
        norm_res = 0.0
      else:
        norm_res = (res.resseq_as_int() - first_resid) / denom
      scaled_index = int(norm_res * len(self._rainbowColors))
      if scaled_index >= len(self._rainbowColors):
        scaled_index = len(self._rainbowColors) - 1
      self._color = self._rainbowColors[scaled_index]

  def get_kin_string(self):
    if self._color is None:
      return ""
    return self._color

  def should_print(self):
    return True


class Range(object):
  """A secondary structure range (HELIX/SHEET/COIL) with sheet linkage info."""

  def __init__(self, ss_type='COIL', chain_id='', init_seq=0, end_seq=0,
               sense=0, strand=1, sheet_id=' ', flipped=False):
    self.type = ss_type
    self.chain_id = chain_id
    self.initSeqNum = init_seq
    self.endSeqNum = end_seq
    self.sense = sense
    self.strand = strand
    self.sheetID = sheet_id
    self.flipped = flipped
    self.previous = None
    self.next = None
    self.duplicateOf = None

  def __str__(self):
    return ('Range: type={}, init={}, end={}, sense={}, strand={}, '
            'sheet={}, flipped={}'.format(
                self.type, self.initSeqNum, self.endSeqNum,
                self.sense, self.strand, self.sheetID, self.flipped))


class RibbonElement(object):
  """A segment of rendered ribbon with start/end spline indices and SS type."""

  def __init__(self, other=None, set_range=None):
    self.start = 0
    self.end = 0
    self.type = 'COIL'
    self.range = set_range
    if other is not None:
      self.like(other)
    if self.range is not None:
      if self.range.type is None or self.range.type == 'TURN':
        self.type = 'COIL'
      else:
        self.type = self.range.type

  def same_sse(self, that):
    return self.type == that.type and self.range == that.range

  def like(self, that):
    self.start = that.start
    self.end = that.end
    self.type = that.type
    self.range = that.range

  def __lt__(self, that):
    a = 0
    b = 0
    if self.range is not None:
      a = self.range.strand
    if that.range is not None:
      b = that.range.strand
    return a < b


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def build_secondary_structure_map(hierarchy, annotation=None, log=None):
  """Build dict[resseq_as_int -> Range] from an annotation object.

  If annotation is None, computes SS via ksdssp. Uses annotation.helices
  and annotation.sheets[].strands directly via get_start_resseq_as_int() /
  get_end_resseq_as_int() instead of parsing PDB-format text records.

  Args:
    hierarchy: iotbx.pdb.hierarchy object
    annotation: iotbx.pdb.secondary_structure.annotation, or None to compute
    log: output stream for messages (may be None)

  Returns:
    dict mapping resseq_as_int -> Range
  """
  import mmtbx.secondary_structure

  if log is None:
    from libtbx.utils import null_out
    log = null_out()

  # Compute SS annotation if not provided
  if annotation is None:
    try:
      params = mmtbx.secondary_structure.manager.get_default_ss_params()
      params.secondary_structure.protein.search_method = "from_ca"
      params = params.secondary_structure
      ss_manager = mmtbx.secondary_structure.manager(
          hierarchy,
          params=params,
          sec_str_from_pdb_file=None,
          log=log)
      annotation = ss_manager.actual_sec_str
    except Exception:
      annotation = None  # fall through to all-COIL

  # Initialize all residues as COIL
  ss_map = {}
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        ss_map[rg.resseq_as_int()] = Range()

  if annotation is None or annotation.is_empty():
    return ss_map

  # Process helices
  for helix in annotation.helices:
    start = helix.get_start_resseq_as_int()
    end = helix.get_end_resseq_as_int()
    if start is None or end is None:
      continue
    r = Range(
        ss_type='HELIX',
        chain_id=helix.start_chain_id.strip() if helix.start_chain_id else '',
        init_seq=start,
        end_seq=end,
        sheet_id=helix.helix_id.strip() if helix.helix_id else '')
    for i in range(start, end + 1):
      ss_map[i] = r

  # Process sheets
  for sheet in annotation.sheets:
    for strand in sheet.strands:
      start = strand.get_start_resseq_as_int()
      end = strand.get_end_resseq_as_int()
      if start is None or end is None:
        continue
      r = Range(
          ss_type='SHEET',
          chain_id=strand.start_chain_id.strip() if strand.start_chain_id else '',
          init_seq=start,
          end_seq=end,
          sense=int(strand.sense) if strand.sense is not None else 0,
          strand=int(strand.strand_id) if strand.strand_id is not None else 1,
          sheet_id=sheet.sheet_id.strip() if sheet.sheet_id else ' ')
      for i in range(start, end + 1):
        ss_map[i] = r

  return ss_map


def consolidate_sheets(secondary_structure):
  """Link sheet strand Range objects by setting previous/next/duplicateOf.

  Args:
    secondary_structure: dict[resseq_as_int -> Range] from
        build_secondary_structure_map()
  """
  unique_strands = {}

  for rng in secondary_structure.values():
    if rng.type != 'SHEET':
      continue

    key = str(int(rng.initSeqNum)) + rng.chain_id
    if key not in unique_strands:
      unique_strands[key] = rng
    else:
      rng.duplicateOf = unique_strands[key]

    for rng2 in secondary_structure.values():
      if rng2.type != 'SHEET':
        continue
      if rng2.sheetID == rng.sheetID and rng2.strand == rng.strand - 1:
        rng.previous = rng2
        rng2.next = rng

  for rng in secondary_structure.values():
    if rng.type != 'SHEET':
      continue
    if rng.duplicateOf is None and rng.next is not None and rng.next.duplicateOf is not None:
      rng.next.duplicateOf.previous = rng
    if rng.duplicateOf is None and rng.previous is not None and rng.previous.duplicateOf is not None:
      rng.previous = rng.previous.duplicateOf


def spline_interpolate(points, n_intervals):
  """Interpolate points using NRUBS B-spline.

  Args:
    points: list of scitbx.matrix.col (3D positions)
    n_intervals: number of intervals between each pair of guide points

  Returns:
    list of scitbx.matrix.col
  """
  nrubs = NRUBS()
  triples = []
  for pt in points:
    triples.append(Triple(pt[0], pt[1], pt[2]))
  result = nrubs.spline(triples, n_intervals)
  return [col((r.x, r.y, r.z)) for r in result]


def get_point_id(start_guide, end_guide, interval, n_intervals,
                 last_point_id, dna_style=False):
  """Generate kinemage point label for a ribbon point.

  Args:
    start_guide: GuidePoint at the start of this segment
    end_guide: GuidePoint at the end of this segment
    interval: current interval index within the segment
    n_intervals: total intervals per segment
    last_point_id: previous point ID string (for deduplication)
    dna_style: if True, use prevRes for early intervals (DNA/RNA)

  Returns:
    tuple of (point_id_str, updated_last_point_id)
  """
  res = start_guide.nextRes
  if dna_style and interval <= n_intervals // 2:
    res = start_guide.prevRes

  buf = (res.atom_groups()[0].resname.strip() + " " +
         res.parent().id.strip() + " " +
         res.resseq.strip() + res.icode)
  point_id = buf.lower().strip()

  if point_id == last_point_id:
    return '"', last_point_id
  else:
    return point_id, point_id


def format_ribbon_point(guides, splines, i, crayon, last_point_id,
                        line_break=False, dna_style=False):
  """Format one spline point as kinemage text.

  Args:
    guides: list of GuidePoint objects
    splines: list of spline-interpolated col points
    i: index into splines
    crayon: Crayon object for coloring
    last_point_id: previous point ID string
    line_break: if True, emit a 'P' (pen-up) marker
    dna_style: if True, use DNA-style point ID labeling

  Returns:
    tuple of (kin_text, updated_last_point_id)
  """
  n_intervals = (len(splines) - 1) // (len(guides) - 3)
  interval = i % n_intervals
  start_idx = (i // n_intervals) + 1

  point_id, last_point_id = get_point_id(
      guides[start_idx], guides[start_idx + 1],
      interval, n_intervals, last_point_id, dna_style)

  ret = "{" + point_id + "}"
  if line_break:
    ret += "P "
  crayon.for_ribbon(guides[start_idx])
  ret += crayon.get_kin_string()
  ret += " "
  pt = splines[i]
  ret += "{:.3f} {:.3f} {:.3f}\n".format(pt[0], pt[1], pt[2])

  return ret, last_point_id


def render_fancy_ribbon(guides, secondary_structure, n_intervals=4,
                        width_alpha=2.0, width_beta=2.2, width_coil=1.0,
                        list_alpha="", list_beta="", list_coil="",
                        list_coil_outline=None, do_rainbow=False,
                        dna_style=False):
  """Render complete ribbon for one contiguous segment as kinemage text.

  This is the core ribbon rendering function, extracted from
  Program.printFancyRibbon() in mmtbx/programs/ribbons.py.

  Args:
    guides: list of GuidePoint objects for this segment
    secondary_structure: dict[resseq_as_int -> Range]
    n_intervals: spline subdivisions per guide segment
    width_alpha: half-width of alpha helix ribbon
    width_beta: half-width of beta sheet ribbon
    width_coil: width of coil
    list_alpha: kinemage list attributes for helix
    list_beta: kinemage list attributes for sheet
    list_coil: kinemage list attributes for coil
    list_coil_outline: kinemage list attributes for coil outline (or None)
    do_rainbow: if True, color by rainbow instead of SS type
    dna_style: if True, use DNA-style point labeling

  Returns:
    kinemage text string
  """
  ret = ""
  sec_struct = secondary_structure

  # Seven strands of guide points: coil, +/-alpha, +/-beta, +/-beta arrowheads
  half_widths = [0.0, -width_alpha / 2, width_alpha / 2,
                 -width_beta / 2, width_beta / 2,
                 -width_beta, width_beta]
  strands = [[] for _ in range(len(half_widths))]
  for g in guides:
    for i in range(len(half_widths)):
      strands[i].append(g.pos + g.dvec * half_widths[i])

  # Seven strands of interpolated points
  splinepts = []
  for i in range(len(strands)):
    splinepts.append(spline_interpolate(strands[i], n_intervals))

  # Discover ribbon elements
  ribbon_elements = []
  rib_element = RibbonElement()
  prev_rib_elt = RibbonElement()
  ribbon_elements.append(rib_element)
  rib_element.type = None

  for i in range(len(guides) - 1 - 2):
    g1 = guides[i + 1]
    g2 = guides[i + 2]

    curr_ss = RibbonElement(set_range=sec_struct.get(g1.nextRes.resseq_as_int(), Range()))
    next_ss = RibbonElement(set_range=sec_struct.get(g2.nextRes.resseq_as_int(), Range()))

    if rib_element.type is None:
      rib_element.like(curr_ss)

    if i == 0 or not rib_element.same_sse(curr_ss):
      if curr_ss.type == 'HELIX' or curr_ss.type == 'SHEET':
        rib_element.end = n_intervals * i + 1
        rib_element = RibbonElement(curr_ss)
        ribbon_elements.append(rib_element)
        rib_element.start = n_intervals * i + 1

    if not rib_element.same_sse(next_ss):
      if curr_ss.type == 'HELIX' or curr_ss.type == 'SHEET':
        end = n_intervals * i + 0
        if curr_ss.type == 'SHEET':
          end += n_intervals - 1
        rib_element.end = end
        rib_element = RibbonElement()
        ribbon_elements.append(rib_element)
        rib_element.type = 'COIL'
        rib_element.start = end

  rib_element.end = len(splinepts[0]) - 1

  # Crayons for coloring
  normal_crayon = Crayon(do_rainbow)
  edge_crayon = Crayon(False)  # No color; deadblack set on vectorlist

  # Sort by strand number
  ribbon_elements.sort()

  # Create list of just sheets for searching
  sheet_elts = [r for r in ribbon_elements if r.type == 'SHEET']

  last_point_id = ""

  for rib_elt in ribbon_elements:
    st_guide = (rib_elt.start // n_intervals) + 2
    end_guide = (rib_elt.end // n_intervals) - 1

    if rib_elt.type == 'HELIX':
      k = (rib_elt.start + rib_elt.end) // 2
      pt = splinepts[2][k]
      v1 = splinepts[1][k] - pt
      v2 = splinepts[1][k + 1] - pt
      cross = v1.cross(v2)
      # Bug fix: original used np.dot(np.linalg.norm(cross), np.linalg.norm(...))
      # which is scalar*scalar, not a dot product. The intent is a sign check
      # for sidedness using the cross product dotted with the guide cvec.
      guide_idx = k // n_intervals + 1
      if guide_idx < len(guides):
        dot = cross.dot(col(guides[guide_idx].cvec))
      else:
        dot = cross.dot(col(guides[-1].cvec))

      crayon = normal_crayon
      ret += "@ribbonlist {fancy helix} " + list_alpha + "\n"
      last_point_id = ""

      for i in range(rib_elt.start, rib_elt.end):
        if dot > 0:
          text, last_point_id = format_ribbon_point(
              guides, splinepts[2], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
          text, last_point_id = format_ribbon_point(
              guides, splinepts[1], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
        else:
          text, last_point_id = format_ribbon_point(
              guides, splinepts[1], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
          text, last_point_id = format_ribbon_point(
              guides, splinepts[2], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
      ret += text

      crayon = edge_crayon
      ret += "@vectorlist {fancy helix edges} width=1 " + list_alpha + " color= deadblack\n"
      last_point_id = ""
      # Left edge
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.start, crayon, last_point_id,
          line_break=True, dna_style=dna_style)
      ret += text
      for i in range(rib_elt.start, rib_elt.end):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[1], i, crayon, last_point_id, dna_style=dna_style)
        ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
      ret += text
      # Right edge
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.start, crayon, last_point_id,
          line_break=True, dna_style=dna_style)
      ret += text
      for i in range(rib_elt.start, rib_elt.end):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[2], i, crayon, last_point_id, dna_style=dna_style)
        ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
      ret += text

    elif rib_elt.type == 'SHEET':
      dot = 0.0

      if rib_elt.range is not None and rib_elt.range.strand != 1:
        for se in sheet_elts:
          if se.range == rib_elt.range.previous:
            prev_rib_elt = se
            break

        prev_range = prev_rib_elt.range
        if prev_range is None:
          prev_rib_elt = rib_elt

        prev_st_guide = (prev_rib_elt.start // n_intervals) + 2
        prev_end_guide = (prev_rib_elt.end // n_intervals) + 1
        cur_closest = st_guide
        prev_closest = prev_st_guide
        close_dist = 1e10
        for i in range(st_guide, end_guide + 2 + 1):
          for j in range(prev_st_guide, prev_end_guide + 2 + 1):
            try:
              O = _FindNamedAtomInResidue(guides[i].prevRes, ["O"])
              N = _FindNamedAtomInResidue(guides[j].prevRes, ["N"])
              if O is None or N is None:
                continue
              # Bug fix: original used O.pos which doesn't exist; use col(O.xyz)
              dist = (col(O.xyz) - col(N.xyz)).length()
              if dist < close_dist:
                close_dist = dist
                cur_closest = i
                prev_closest = j
            except Exception:
              pass

        k_cur = min(n_intervals * (cur_closest - 1), len(splinepts[4]) - 1)
        k_prev = min(n_intervals * (prev_closest - 1), len(splinepts[4]) - 1)
        pt_cur = splinepts[4][k_cur]
        v1_cur = splinepts[3][k_cur] - pt_cur
        if k_cur + 1 < len(splinepts[3]):
          v2_cur = splinepts[3][k_cur + 1] - pt_cur
          cross_cur = v1_cur.cross(v2_cur)
          pt_prev = splinepts[4][k_prev]
          v1_prev = splinepts[3][k_prev] - pt_prev
          if k_prev + 1 < len(splinepts[3]):
            v2_prev = splinepts[3][k_prev + 1] - pt_prev
            cross_prev = v1_prev.cross(v2_prev)
            dot = cross_cur.dot(cross_prev)

      crayon = normal_crayon
      ret += "@ribbonlist {fancy sheet} " + list_beta + "\n"
      last_point_id = ""
      for i in range(rib_elt.start, rib_elt.end - 1):
        if ((dot < 0 and not prev_rib_elt.range.flipped) or
            (dot > 0 and prev_rib_elt.range.flipped)):
          text, last_point_id = format_ribbon_point(
              guides, splinepts[4], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
          text, last_point_id = format_ribbon_point(
              guides, splinepts[3], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
          if rib_elt.range is not None:
            rib_elt.range.flipped = True
        else:
          text, last_point_id = format_ribbon_point(
              guides, splinepts[3], i, crayon, last_point_id, dna_style=dna_style)
          ret += text
          text, last_point_id = format_ribbon_point(
              guides, splinepts[4], i, crayon, last_point_id, dna_style=dna_style)
          ret += text

      # Arrowhead
      if ((dot < 0 and not prev_rib_elt.range.flipped) or
          (dot > 0 and prev_rib_elt.range.flipped)):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[6], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
        ret += text
        text, last_point_id = format_ribbon_point(
            guides, splinepts[5], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
        ret += text
        text, last_point_id = format_ribbon_point(
            guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
        ret += text
      else:
        text, last_point_id = format_ribbon_point(
            guides, splinepts[5], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
        ret += text
        text, last_point_id = format_ribbon_point(
            guides, splinepts[6], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
        ret += text
        text, last_point_id = format_ribbon_point(
            guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
        ret += text

      # Borders
      crayon = edge_crayon
      ret += "@vectorlist {fancy sheet edges} width=1 " + list_beta + " color= deadblack\n"
      last_point_id = ""
      # Left edge
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.start, crayon, last_point_id,
          line_break=True, dna_style=dna_style)
      ret += text
      for i in range(rib_elt.start, rib_elt.end - 1):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[3], i, crayon, last_point_id, dna_style=dna_style)
        ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[5], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
      ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
      ret += text
      # Right edge
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.start, crayon, last_point_id,
          line_break=True, dna_style=dna_style)
      ret += text
      for i in range(rib_elt.start, rib_elt.end - 1):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[4], i, crayon, last_point_id, dna_style=dna_style)
        ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[6], rib_elt.end - 2, crayon, last_point_id, dna_style=dna_style)
      ret += text
      text, last_point_id = format_ribbon_point(
          guides, splinepts[0], rib_elt.end, crayon, last_point_id, dna_style=dna_style)
      ret += text

    else:  # Coil
      if list_coil_outline is not None:
        crayon = normal_crayon
        ret += "@vectorlist {fancy coil edges} " + list_coil_outline + "\n"
        last_point_id = ""
        for i in range(rib_elt.start, rib_elt.end + 1):
          text, last_point_id = format_ribbon_point(
              guides, splinepts[0], i, crayon, last_point_id, dna_style=dna_style)
          ret += text

      crayon = normal_crayon
      ret += "@vectorlist {fancy coil} " + list_coil + "\n"
      last_point_id = ""
      for i in range(rib_elt.start, rib_elt.end + 1):
        text, last_point_id = format_ribbon_point(
            guides, splinepts[0], i, crayon, last_point_id, dna_style=dna_style)
        ret += text

    crayon = normal_crayon

  return ret


def generate_chain_ribbons(chain, secondary_structure, chain_id,
                           chain_color="white", do_protein=True,
                           do_nucleic_acid=True, untwist=True,
                           dna_style=False, do_plain_coils=False,
                           coil_width=1.0, color_by="secondary_structure",
                           nucleic_acid_as_helix=True,
                           has_dna=False, has_rna=False):
  """Generate complete ribbon kinemage text for one chain.

  Top-level function called per chain to produce ribbon content including
  colorset definitions, subgroup headers, and all ribbon lists.

  Args:
    chain: iotbx.pdb.hierarchy chain object
    secondary_structure: dict[resseq_as_int -> Range]
    chain_id: chain identifier string
    chain_color: base color for this chain
    do_protein: generate protein ribbons
    do_nucleic_acid: generate nucleic acid ribbons
    untwist: remove excess twist from ribbons
    dna_style: use DNA-style ribbon orientation
    do_plain_coils: if True, omit coil outlines (halo)
    coil_width: width of the coil part
    color_by: "rainbow", "secondary_structure", or "solid"
    nucleic_acid_as_helix: treat nucleic acids as helix in SS map
    has_dna: whether DNA is present in the model
    has_rna: whether RNA is present in the model

  Returns:
    kinemage text string for this chain's ribbons
  """
  out = ""
  do_rainbow = (color_by == "rainbow")

  # --- Protein ribbons ---
  if do_protein:
    contiguous_lists = find_contiguous_protein_residues(chain)
    if len(contiguous_lists) > 0:
      out += "@subgroup {ribbon %s}\n" % chain_id

      if chain_color == "white" and color_by != "solid":
        out += "@colorset {{alph{}}} red\n".format(chain_id)
        out += "@colorset {{beta{}}} lime\n".format(chain_id)
        out += "@colorset {{coil{}}} white\n".format(chain_id)
      else:
        out += "@colorset {{alph{}}} {}\n".format(chain_id, chain_color)
        out += "@colorset {{beta{}}} {}\n".format(chain_id, chain_color)
        out += "@colorset {{coil{}}} {}\n".format(chain_id, chain_color)

      for contig in contiguous_lists:
        guidepoints = make_protein_guidepoints(contig)
        if untwist:
          untwist_ribbon(guidepoints)
        if do_plain_coils:
          out += render_fancy_ribbon(
              guidepoints, secondary_structure,
              width_alpha=2.0, width_beta=2.2, width_coil=coil_width,
              list_alpha="color= {alph" + chain_id + "} master= {protein} master= {ribbon} master= {alpha}",
              list_beta="color= {beta" + chain_id + "} master= {protein} master= {ribbon} master= {beta}",
              list_coil="width= 4 color= {coil" + chain_id + "} master= {protein} master= {ribbon} master= {coil}",
              do_rainbow=do_rainbow, dna_style=dna_style)
        else:
          out += render_fancy_ribbon(
              guidepoints, secondary_structure,
              width_alpha=2.0, width_beta=2.2, width_coil=coil_width,
              list_alpha="color= {alph" + chain_id + "} master= {protein} master= {ribbon} master= {alpha}",
              list_beta="color= {beta" + chain_id + "} master= {protein} master= {ribbon} master= {beta}",
              list_coil="width= 4 fore color= {coil" + chain_id + "} master= {protein} master= {ribbon} master= {coil}",
              list_coil_outline="width= 6 rear color= deadblack master= {protein} master= {ribbon} master= {coil}",
              do_rainbow=do_rainbow, dna_style=dna_style)

  # --- Nucleic acid ribbons ---
  if do_nucleic_acid:
    contiguous_lists = find_contiguous_nucleic_acid_residues(chain)
    if len(contiguous_lists) > 0:
      out += "@subgroup {ribbon %s NA}\n" % chain_id

      if chain_color == "white":
        out += "@colorset {{nucl{}}} lime\n".format(chain_id)
        out += "@colorset {{ncoi{}}} white\n".format(chain_id)
      else:
        out += "@colorset {{nucl{}}} {}\n".format(chain_id, chain_color)
        out += "@colorset {{ncoi{}}} {}\n".format(chain_id, chain_color)

      for contig in contiguous_lists:
        guidepoints = make_nucleic_acid_guidepoints(contig)
        if untwist:
          untwist_ribbon(guidepoints)
        use_dna_style = dna_style or (has_dna and has_rna and chain_has_DNA(chain))
        if use_dna_style:
          swap_edge_and_face(guidepoints)

        out += render_fancy_ribbon(
            guidepoints, secondary_structure,
            width_alpha=3.0, width_beta=3.0, width_coil=coil_width,
            list_alpha="color= {nucl" + chain_id + "} master= {nucleic acid} master= {ribbon} master= {RNA helix?}",
            list_beta="color= {nucl" + chain_id + "} master= {nucleic acid} master= {ribbon} master= {A-form}",
            list_coil="width= 4 color= {ncoi" + chain_id + "} master= {nucleic acid} master= {ribbon} master= {coil}",
            do_rainbow=do_rainbow, dna_style=use_dna_style)

  return out
