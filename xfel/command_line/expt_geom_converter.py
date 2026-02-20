# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.expt_geom_converter
#!/usr/bin/env python3
"""
Convert between DIALS .expt files and CrystFEL .geom files.

Coordinate conventions
----------------------
dxtbx / imgCIF:  +x right, +y up, +z toward source (opposite beam).
CrystFEL:        +x right, +y up, +z along beam (opposite dxtbx z).

Conversion (x,y unchanged; z sign-flipped):
  corner_x = origin_x_mm / pixel_size_mm      (pixels from beam axis)
  corner_y = origin_y_mm / pixel_size_mm
  z_m      = -origin_z_mm / 1000               (metres, CrystFEL sign)
  coffset  = z_m - clen_m
"""

import sys
import math
import os
import json
import re
import tempfile
from types import SimpleNamespace
from libtbx.phil import parse as phil_parse

from dxtbx.model.experiment_list import ExperimentList


phil_scope = phil_parse("""
mode = *expt2geom geom2expt validate
 .type = choice
 .help = Operation to perform.

expt2geom {
 input_expt = None
  .type = path
  .help = Input DIALS experiment file (.expt).
 output_geom = None
  .type = path
  .help = Output CrystFEL geometry file (.geom).
 clen = None
  .type = float
  .help = Override nominal camera length in metres (default: mean panel z).
 data_path = /entry/data/data
  .type = str
  .help = HDF5 path to image data written into the .geom file.
 adu_per_ev = None
  .type = float
  .help = ADU per eV. Auto-computed from detector gain if not set.
 max_adu = None
  .type = int
  .help = Maximum trusted ADU. Auto-computed from trusted_range if not set.
 mask_path = None
  .type = str
  .help = HDF5 path to pixel mask (optional).
}

geom2expt {
 input_geom = None
  .type = path
  .help = Input CrystFEL geometry file (.geom).
 template_expt = None
  .type = path
  .help = Template DIALS experiment file (.expt).
 output_expt = None
  .type = path
  .help = Output DIALS experiment file (.expt).
 wavelength = None
  .type = float
  .help = Override beam wavelength in Angstroms.
}

validate {
 input_expt = None
  .type = path
  .help = Input DIALS experiment file (.expt) to round-trip validate.
 tol_origin = 0.001
  .type = float
  .help = Origin displacement tolerance in mm.
 tol_axis = 0.001
  .type = float
  .help = Axis angle tolerance in degrees.
}
""")


# --- Vector helpers ---

def format_vec2(fx, fy, threshold=1e-6):
  """Format a CrystFEL 2D direction vector, e.g. '0.999000x +0.001000y'."""
  parts = []
  if abs(fx) > threshold:
    parts.append(f"{fx:.6f}x")
  if abs(fy) > threshold:
    # Include '+' sign only when there is already an x term
    if parts:
      parts.append(f"{fy:+.6f}y")
    else:
      parts.append(f"{fy:.6f}y")
  if not parts:
    parts = ["0.000000x"]
  return " ".join(parts)


def _dot(a, b):
  """Dot product of two equal-length sequences."""
  return sum(x * y for x, y in zip(a, b))


def _angle_deg_2d(a, b):
  """Angle in degrees between two vectors, using only their x,y components."""
  a2 = (a[0], a[1])
  b2 = (b[0], b[1])
  na = math.sqrt(_dot(a2, a2))
  nb = math.sqrt(_dot(b2, b2))
  if na < 1e-12 or nb < 1e-12:
    return 0.0
  cos_angle = _dot(a2, b2) / (na * nb)
  cos_angle = max(-1.0, min(1.0, cos_angle))
  return math.degrees(math.acos(cos_angle))


# --- expt -> geom ---

def load_expt(expt_file):
  """Load experiment 0 from a .expt file.  Returns (beam, detector)."""
  expts = ExperimentList.from_file(expt_file, check_format=False)
  if len(expts) == 0:
    sys.exit(f"ERROR: no experiments found in {expt_file}")

  return expts[0].beam, expts[0].detector


def extract_panel_data(detector, beam, clen_override=None):
  """Extract per-panel CrystFEL geometry from a dxtbx detector.

  Returns (panels_data, clen_m) where clen_m is in metres and each panel
  dict has keys: name, min/max_fs/ss, corner_x/y, coffset, fs_x/y,
  ss_x/y, res, z_m.
  """
  # Warn if beam is significantly off-axis.
  unit_s0 = beam.get_unit_s0()
  beam_dir = (-unit_s0[0], -unit_s0[1], -unit_s0[2])
  tilt = math.sqrt(beam_dir[0] ** 2 + beam_dir[1] ** 2)
  if tilt > 0.01:
    print(
      f"WARNING: beam is significantly off-axis "
      f"(unit_s0 x={unit_s0[0]:.4f}, y={unit_s0[1]:.4f}). "
      "corner_x/y assume beam on axis; geometry may be inaccurate."
    )

  panels_data = []

  for i, panel in enumerate(detector):
    origin     = panel.get_origin()
    fast       = panel.get_fast_axis()
    slow       = panel.get_slow_axis()
    pixel_size = panel.get_pixel_size()
    image_size = panel.get_image_size()
    name       = panel.get_name() or f"p{i}"

    # CrystFEL reserves 'bad*' and 'rigid_group*' prefixes.
    if name.lower().startswith(("bad", "rigid_group")):
      name = f"p{i}"

    pixel_size_mm = pixel_size[0]

    if abs(pixel_size[0] - pixel_size[1]) > 1e-9:
      print(
        f"WARNING: panel {name} has non-square pixels "
        f"({pixel_size[0]:.6f} x {pixel_size[1]:.6f} mm). "
        "CrystFEL uses a single res value; using fast-axis pixel size."
      )

    res      = 1.0 / (pixel_size_mm * 1e-3)   # pixels per metre
    corner_x = origin[0] / pixel_size_mm       # pixels from beam axis
    corner_y = origin[1] / pixel_size_mm
    z_m      = -origin[2] / 1000.0             # sign flip, dxtbx -> CrystFEL

    fs_x, fs_y, fs_z = fast
    ss_x, ss_y, ss_z = slow

    if abs(fs_z) > 0.01 or abs(ss_z) > 0.01:
      print(
        f"WARNING: panel {name} has non-negligible z components in "
        f"fast axis (z={fs_z:.4f}) or slow axis (z={ss_z:.4f}). "
        "CrystFEL requires panels perpendicular to the beam; "
        "z components will be dropped and geometry may be approximate."
      )

    panels_data.append(dict(
      name=name,
      min_fs=0,
      max_fs=image_size[0] - 1,
      min_ss=0,
      max_ss=image_size[1] - 1,
      corner_x=corner_x,
      corner_y=corner_y,
      coffset=0.0,  # computed after clen is known
      fs_x=float(fs_x),
      fs_y=float(fs_y),
      ss_x=float(ss_x),
      ss_y=float(ss_y),
      res=res,
      z_m=z_m,
    ))

  # Nominal camera length
  if clen_override is not None:
    clen_m = clen_override
    print(f"Using user-specified clen = {clen_m:.6f} m")
  else:
    z_values = [pd["z_m"] for pd in panels_data]
    clen_m = sum(z_values) / len(z_values)
    print(
      f"Computed clen = {clen_m * 1000:.3f} mm "
      f"(mean z of {len(panels_data)} panel(s))"
    )

  # Compute coffset for each panel
  for pd in panels_data:
    pd["coffset"] = pd["z_m"] - clen_m

  return panels_data, clen_m


def _auto_adu_per_ev(detector, wavelength_A):
  """Compute adu_per_ev from detector gain and wavelength.
  Returns the value, or None if gain is unavailable."""
  gains = [panel.get_gain() for panel in detector]
  if not gains or gains[0] is None:
    return None
  mean_gain = sum(gains) / len(gains)
  if max(gains) != min(gains):
    print(f"WARNING: non-uniform gain across panels; using mean {mean_gain:.6f}")
  photon_energy_eV = 12398.4 / wavelength_A
  return mean_gain / photon_energy_eV


def _auto_max_adu(detector):
  """Compute max_adu from the upper bound of panel trusted_range.
  Returns the value as int, or None if unavailable."""
  uppers = []
  for panel in detector:
    tr = panel.get_trusted_range()
    if tr is not None and len(tr) >= 2:
      uppers.append(tr[1])
  if not uppers:
    return None
  result = min(uppers)  # conservative: use lowest upper bound
  if max(uppers) != min(uppers):
    print(f"WARNING: non-uniform trusted_range across panels; using minimum {int(result)}")
  return int(result)


def write_geom(panels_data, clen_m, wavelength_A, geom_file, args):
  """Write a CrystFEL .geom file from extracted panel data."""
  photon_energy_eV = 12398.4 / wavelength_A
  n_panels = len(panels_data)
  panel_names = [pd["name"] for pd in panels_data]

  # All panels must share the same res for the single global line.
  res_global = panels_data[0]["res"]
  if any(abs(pd["res"] - res_global) > 1.0 for pd in panels_data):
    print(
      "WARNING: panels have different res values; "
      f"writing global res = {res_global:.2f} (from first panel). "
      "Review the output .geom file."
    )

  with open(geom_file, "w") as f:
    p = lambda *a, **kw: print(*a, file=f, **kw)

    # Header
    p("; CrystFEL geometry file")
    p("; Generated by expt_geom_converter.py")
    p(";")
    p("; IMPORTANT: review before use:")
    p(";   1. photon_energy / adu_per_eV")
    p(";   2. data path and dim entries (match your HDF5 layout)")
    p(";   3. mask path if applicable")
    p(";   4. clen (nominal detector distance in metres)")
    p(";")
    p()

    # Global parameters
    p("; ===== Global parameters =====")
    p()
    p(f"photon_energy = {photon_energy_eV:.2f}")
    p(f";   (wavelength = {wavelength_A:.6f} Ang)")
    p()

    if args.adu_per_ev:
      p(f"adu_per_eV = {args.adu_per_ev:.8f}")
      p()
    else:
      p(f"; adu_per_eV = ???  ; SET THIS: detector ADU per eV")
      p(f"; Example for 1 ADU/photon: adu_per_eV = "
       f"{1.0 / photon_energy_eV:.8f}")
      p()

    if args.max_adu:
      p(f"max_adu = {args.max_adu}")
      p()

    p(f"clen = {clen_m:.6f}")
    p(f";   ({clen_m * 1000:.3f} mm)")
    p()

    # Global res -- written once, not per-panel.
    p(f"res = {res_global:.2f}")
    p()

    # Data layout
    p("; ===== Data layout =====")
    p("; Adjust these to match your HDF5 file structure.")
    p()
    p(f"data = {args.data_path}")
    p("dim0 = %")
    if n_panels > 1:
      p("dim1 = %")
      p("dim2 = ss")
      p("dim3 = fs")
    else:
      p("dim1 = ss")
      p("dim2 = fs")
    p()

    # Mask
    if args.mask_path:
      p("; ===== Mask =====")
      p()
      p(f"mask = {args.mask_path}")
      p("mask_good = 0x0")
      p("mask_bad  = 0x1")
      p()

    # Rigid groups (one per panel)
    p("; ===== Rigid groups =====")
    p("; One rigid group per panel.")
    for i, name in enumerate(panel_names):
      p(f"rigid_group_m{i} = {name}")
    p()
    all_groups = ", ".join(f"m{i}" for i in range(n_panels))
    p(f"rigid_group_collection_modules = {all_groups}")
    p()

    # Panel definitions
    p("; ===== Panel definitions =====")
    p()
    for i, pd in enumerate(panels_data):
      name = pd["name"]
      p(f"; --- Panel {i}: {name} ---")
      p(f"{name}/min_fs = {pd['min_fs']}")
      p(f"{name}/max_fs = {pd['max_fs']}")
      p(f"{name}/min_ss = {pd['min_ss']}")
      p(f"{name}/max_ss = {pd['max_ss']}")
      p(f"{name}/fs = {format_vec2(pd['fs_x'], pd['fs_y'])}")
      p(f"{name}/ss = {format_vec2(pd['ss_x'], pd['ss_y'])}")
      p(f"{name}/corner_x = {pd['corner_x']:.4f}")
      p(f"{name}/corner_y = {pd['corner_y']:.4f}")
      p(f"{name}/coffset = {pd['coffset']:.6f}")
      if n_panels > 1:
        p(f"{name}/dim1 = {i}")
      p()

  print(f"Wrote {n_panels} panel(s) to {geom_file}.")


# --- geom -> expt ---

def parse_vec2(s):
  """Parse a CrystFEL 2D vector string like '0.999000x +0.001000y'."""
  x_val = 0.0
  y_val = 0.0
  for token in s.split():
    token = token.strip()
    if token.endswith("x"):
      x_val = float(token[:-1])
    elif token.endswith("y"):
      y_val = float(token[:-1])
  return x_val, y_val


def _geom_panel_value(flat, panel_name, suffix, default=None):
  """Look up 'panel_name/suffix' in the flat key-value dict from a .geom file."""
  val = flat.get(f"{panel_name}/{suffix}")
  if val is not None:
    return val
  if default is not None:
    return default
  raise ValueError(
    f"Required key '{panel_name}/{suffix}' not found in .geom file"
  )


def parse_geom(geom_file):
  """Parse a .geom file and return panel data in dxtbx coordinates.

  Returns (panels_data, clen_m, wavelength_A).  Each panel dict has keys:
  name, min/max_fs/ss, origin_x/y_mm, origin_z_dxtbx_mm, fast_axis,
  slow_axis.
  """
  # Build a flat key=value dict, stripping comments.
  flat = {}
  with open(geom_file) as fh:
    for line in fh:
      line = line.split(";")[0].strip()
      if "=" not in line:
        continue
      key, _, val = line.partition("=")
      flat[key.strip()] = val.strip()

  if "clen" not in flat:
    raise ValueError("Required global key 'clen' not found in .geom file")
  clen_m = float(flat["clen"])

  if "res" not in flat:
    raise ValueError("Required global key 'res' not found in .geom file")
  pixel_size_mm = 1000.0 / float(flat["res"])

  photon_energy = float(flat.get("photon_energy", "0") or "0")
  wavelength_A = 12398.4 / photon_energy if photon_energy > 0 else 0.0

  # Discover panel names from 'panelname/corner_x' keys.
  panel_names = sorted({
    m.group(1)
    for key in flat
    if (m := re.match(r"^(.+)/corner_x$", key))
  })

  panels_data = []
  for pname in panel_names:
    get = lambda suffix, default=None, _pn=pname: (
      _geom_panel_value(flat, _pn, suffix, default)
    )

    corner_x  = float(get("corner_x"))
    corner_y  = float(get("corner_y"))
    coffset_m = float(get("coffset", "0"))

    fs_x, fs_y = parse_vec2(get("fs"))
    ss_x, ss_y = parse_vec2(get("ss"))

    # Convert to dxtbx coordinates (x,y unchanged; z sign-flipped).
    origin_x_mm           = corner_x * pixel_size_mm
    origin_y_mm           = corner_y * pixel_size_mm
    origin_z_dxtbx_mm     = -(clen_m + coffset_m) * 1000.0

    # Gram-Schmidt: re-orthogonalise slow against fast.  CrystFEL rounds
    # fs/ss independently, which can violate dxtbx's orthogonality check.
    dot_fs = fs_x * ss_x + fs_y * ss_y
    norm2  = fs_x * fs_x + fs_y * fs_y
    if norm2 > 0:
      ss_x -= dot_fs / norm2 * fs_x
      ss_y -= dot_fs / norm2 * fs_y

    panels_data.append(dict(
      name=pname,
      min_fs=int(get("min_fs")),
      max_fs=int(get("max_fs")),
      min_ss=int(get("min_ss")),
      max_ss=int(get("max_ss")),
      origin_x_mm=origin_x_mm,
      origin_y_mm=origin_y_mm,
      origin_z_dxtbx_mm=origin_z_dxtbx_mm,
      fast_axis=(fs_x, fs_y, 0.0),
      slow_axis=(ss_x, ss_y, 0.0),
    ))

  return panels_data, clen_m, wavelength_A


def _reset_hierarchy_to_identity(hierarchy):
  """Set all group-level transforms in a detector hierarchy to identity.

  Leaf nodes (those with a "panel" key) are left unchanged.
  """
  hierarchy["origin"]    = [0, 0, 0]
  hierarchy["fast_axis"] = [1, 0, 0]
  hierarchy["slow_axis"] = [0, 1, 0]
  for child in hierarchy.get("children", []):
    if "panel" not in child:
      child["origin"]    = [0, 0, 0]
      child["fast_axis"] = [1, 0, 0]
      child["slow_axis"] = [0, 1, 0]


def build_expt_from_template(template_expt, panels_data, wavelength_A, output_expt):
  """Replace geometry in a template .expt with parsed .geom panel data.

  Matches panels by name when possible, falls back to index order.
  Resets detector hierarchy to identity so dxtbx uses panel coordinates
  directly.
  """
  with open(template_expt) as fh:
    data = json.load(fh)

  panels = data["detector"][0]["panels"]

  if len(panels) != len(panels_data):
    raise ValueError(
      f"Panel count mismatch: template has {len(panels)} panel(s), "
      f".geom file has {len(panels_data)} panel(s)."
    )

  template_name_to_idx = {p["name"]: j for j, p in enumerate(panels)}
  use_name_match = all(
    pd["name"] in template_name_to_idx for pd in panels_data
  )

  if not use_name_match:
    print(
      "WARNING: panel names in .geom file do not all match template panel "
      "names. Falling back to index-order matching."
    )

  for i, pd in enumerate(panels_data):
    j = template_name_to_idx[pd["name"]] if use_name_match else i
    panels[j]["origin"] = [
      pd["origin_x_mm"],
      pd["origin_y_mm"],
      pd["origin_z_dxtbx_mm"],
    ]
    panels[j]["fast_axis"] = list(pd["fast_axis"])
    panels[j]["slow_axis"] = list(pd["slow_axis"])

  hierarchy = data["detector"][0].get("hierarchy", {})
  _reset_hierarchy_to_identity(hierarchy)
  data["detector"][0]["hierarchy"] = hierarchy

  if wavelength_A > 0:
    data["beam"][0]["wavelength"] = wavelength_A

  with open(output_expt, "w") as fh:
    json.dump(data, fh, indent=2)

  print(f"Wrote updated experiment to {output_expt}.")


# --- CLI ---

def _cmd_expt2geom(params):
  """Entry point for the expt2geom subcommand."""
  input_expt = str(params.expt2geom.input_expt)
  output_geom = str(params.expt2geom.output_geom)
  beam, detector = load_expt(input_expt)
  wavelength_A = beam.get_wavelength()
  print(f"Loaded experiment from {input_expt}  "
     f"({len(detector)} panel(s), λ = {wavelength_A:.6f} Å)")

  if params.expt2geom.adu_per_ev is None:
    auto = _auto_adu_per_ev(detector, wavelength_A)
    if auto is not None:
      print(f"Auto-computed adu_per_ev = {auto:.8f} from detector gain")
      params.expt2geom.adu_per_ev = auto

  if params.expt2geom.max_adu is None:
    auto = _auto_max_adu(detector)
    if auto is not None:
      print(f"Auto-computed max_adu = {auto} from trusted_range")
      params.expt2geom.max_adu = auto

  panels_data, clen_m = extract_panel_data(
    detector, beam, clen_override=params.expt2geom.clen
  )

  write_geom(panels_data, clen_m, wavelength_A, output_geom, params.expt2geom)
  print("Done.")


def _cmd_geom2expt(params):
  """Entry point for the geom2expt subcommand."""
  input_geom = str(params.geom2expt.input_geom)
  template_expt = str(params.geom2expt.template_expt)
  output_expt = str(params.geom2expt.output_expt)
  panels_data, clen_m, wavelength_A = parse_geom(input_geom)
  print(
    f"Parsed {len(panels_data)} panel(s) from {input_geom}  "
    f"(clen = {clen_m * 1000:.3f} mm, λ = {wavelength_A:.6f} Å)"
  )

  if params.geom2expt.wavelength is not None:
    wavelength_A = params.geom2expt.wavelength
    print(f"Using user-specified wavelength = {wavelength_A:.6f} Å")

  build_expt_from_template(
    template_expt, panels_data, wavelength_A, output_expt
  )
  print("Done.")


def validate(expt_file, tol_origin_mm=0.001, tol_axis_deg=0.001):
  """Round-trip validation: expt -> geom -> expt -> compare.

  Returns True if all panels match within tolerances.
  """
  write_stub = SimpleNamespace(
    adu_per_ev=None,
    max_adu=None,
    mask_path=None,
    data_path="/entry/data/data",
  )

  with tempfile.TemporaryDirectory() as tmpdir:
    tmp_geom = tmpdir + "/tmp.geom"
    tmp_expt = tmpdir + "/tmp.expt"

    # Forward pass: expt -> geom
    beam, det_orig = load_expt(expt_file)
    panels_data, clen_m = extract_panel_data(det_orig, beam)
    write_geom(panels_data, clen_m, beam.get_wavelength(), tmp_geom, write_stub)

    # Reverse pass: geom -> expt
    rt_panels, rt_clen, rt_wl = parse_geom(tmp_geom)
    build_expt_from_template(expt_file, rt_panels, rt_wl, tmp_expt)

    # Load round-tripped detector
    _, det_rt = load_expt(tmp_expt)

  # Compare panel-by-panel
  col_w = [22, 13, 13, 13, 8]
  header = (
    f"{'Panel':<{col_w[0]}}"
    f"{'Origin(mm)':<{col_w[1]}}"
    f"{'Fast(deg)':<{col_w[2]}}"
    f"{'Slow(deg)':<{col_w[3]}}"
    f"{'Status':<{col_w[4]}}"
  )
  print(header)
  print("-" * sum(col_w))

  results = []
  for panel_orig, panel_rt in zip(det_orig, det_rt):
    name = panel_orig.get_name()

    orig_o = panel_orig.get_origin()
    rt_o   = panel_rt.get_origin()
    d_origin = math.sqrt(sum((a - b) ** 2 for a, b in zip(orig_o, rt_o)))

    orig_fast = panel_orig.get_fast_axis()
    rt_fast   = panel_rt.get_fast_axis()
    d_fast = _angle_deg_2d(orig_fast, rt_fast)

    orig_slow = panel_orig.get_slow_axis()
    rt_slow   = panel_rt.get_slow_axis()
    d_slow = _angle_deg_2d(orig_slow, rt_slow)

    passed = (
      d_origin <= tol_origin_mm
      and d_fast <= tol_axis_deg
      and d_slow <= tol_axis_deg
    )
    status = "PASS" if passed else "FAIL"
    results.append((name, d_origin, d_fast, d_slow, passed))

    print(
      f"{name:<{col_w[0]}}"
      f"{d_origin:<{col_w[1]}.6f}"
      f"{d_fast:<{col_w[2]}.6f}"
      f"{d_slow:<{col_w[3]}.6f}"
      f"{status:<{col_w[4]}}"
    )

  n_pass = sum(1 for r in results if r[4])
  n_total = len(results)
  max_origin = max(r[1] for r in results) if results else 0.0
  max_axis   = max(max(r[2], r[3]) for r in results) if results else 0.0

  # Check for z-component loss (inherent to .geom format)
  max_z_fast = max(abs(det_orig[i].get_fast_axis()[2]) for i in range(len(det_orig)))
  max_z_slow = max(abs(det_orig[i].get_slow_axis()[2]) for i in range(len(det_orig)))

  print()
  print(f"Summary: {n_pass}/{n_total} panels PASS")
  print(f"  Max origin error : {max_origin:.6f} mm  (tol={tol_origin_mm} mm)")
  print(f"  Max axis error   : {max_axis:.6f} deg (tol={tol_axis_deg} deg, x-y plane only)")
  if max_z_fast > 1e-4 or max_z_slow > 1e-4:
    print(f"  Note: axis z-components lost in .geom round-trip "
       f"(max |z|: fast={max_z_fast:.6f}, slow={max_z_slow:.6f})")

  return n_pass == n_total


def _cmd_validate(params):
  """Entry point for the validate subcommand."""
  passed = validate(
    params.validate.input_expt,
    params.validate.tol_origin,
    params.validate.tol_axis,
  )
  sys.exit(0 if passed else 1)


def parse_phil_args(argv):
  """Parse command-line arguments as phil expressions or phil files."""
  if "-h" in argv or "--help" in argv:
    print("Usage: python expt_geom_converter.py [params.phil] [key=value ...]")
    print()
    print("Modes: expt2geom, geom2expt, validate")
    print()
    print("Use -c to show full parameter definitions.")
    sys.exit(0)
  if "-c" in argv:
    phil_scope.show(attributes_level=2)
    sys.exit(0)

  user_phil = []
  for arg in argv:
    if os.path.isfile(arg):
      user_phil.append(phil_parse(file_name=arg))
    else:
      try:
        user_phil.append(phil_parse(arg))
      except Exception:
        sys.exit(f"ERROR: unrecognized argument: {arg}")

  return phil_scope.fetch(sources=user_phil).extract()


def main():
  params = parse_phil_args(sys.argv[1:])

  if params.mode == "expt2geom":
    _cmd_expt2geom(params)
  elif params.mode == "geom2expt":
    _cmd_geom2expt(params)
  elif params.mode == "validate":
    _cmd_validate(params)


if __name__ == "__main__":
  main()
