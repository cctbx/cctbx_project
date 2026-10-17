""" Convert a PDB or mmCIF file into a SHELX .ins that preserves residues.

  iotbx.python pdb_as_shelx_ins.py <in.pdb|in.cif> [out.ins] [--p1] [--no-water]

Why not `iotbx.shelx.writer.generator`: it refuses duplicate atom names, and a
protein has one `CA` per residue; it emits no RESI records; and it sorts
scatterers by element, which destroys residue grouping. Residue grouping is the
whole point here, because Olex2 carries residues as SHELX `RESI` blocks
(`xlib/residue.h`, `TAsymmUnit::ResidueRegistry`) and fragHAR is residue-based.
The cell, LATT and SFAC emission below therefore overlaps that module by a
dozen lines; the alternative was to make `generator` residue-aware, which would
change a function other code already depends on.

Why not rely on Olex2's own PDB reader: `xlib/pdb.cpp` handles only CRYST1,
ATOM and ANISOU. It ignores HETATM, so ligands, metals and waters vanish;
it ignores TER, so chain breaks are invisible; and it allocates but never reads
the insertion code. Converting once, deliberately, loses less.

The output is a normal refinable model: bare coordinates and Uiso, occupancy
fixed at 1 in the SHELX manner (11.00000), one RESI block per residue.
"""
from __future__ import absolute_import, division, print_function

import math
import os
import sys

# SHELX puts C and H first by convention, then ascending Z
PRIOR_ELEMENTS = ("C", "H")


def sanitise_label(name, used):
  """A SHELX label: at most four characters, no spaces or primes.

  Uniqueness is only required within a residue, which is what RESI scoping
  buys, so a collision here is rare and is resolved by numbering rather than by
  renaming the atom to something unrecognisable.
  """
  s = name.strip().replace("'", "").replace('"', "").replace(" ", "")
  if not s:
    s = "X"
  s = s[:4]
  if s not in used:
    used.add(s)
    return s
  for i in range(1, 100):
    c = ("%s%d" % (s[:3], i))[:4]
    if c not in used:
      used.add(c)
      return c
  used.add(s)
  return s


def box_symmetry(hierarchy, margin=10.0):
  """A P1 box around the model, for entries with no CRYST1 (NMR, predictions)."""
  from cctbx import crystal
  xyz = hierarchy.atoms().extract_xyz()
  lo = [min(p[i] for p in xyz) for i in range(3)]
  hi = [max(p[i] for p in xyz) for i in range(3)]
  dims = [hi[i] - lo[i] + 2*margin for i in range(3)]
  return crystal.symmetry(
    unit_cell=(dims[0], dims[1], dims[2], 90, 90, 90),
    space_group_symbol="P 1"), lo, margin


def convert(in_path, out_path, force_p1=False, keep_water=True):
  import iotbx.pdb
  from cctbx import sgtbx
  from cctbx.eltbx import tiny_pse

  inp = iotbx.pdb.input(file_name=in_path)
  hierarchy = inp.construct_hierarchy()
  hierarchy.remove_alt_confs(always_keep_one_conformer=True)

  cs = inp.crystal_symmetry()
  offset = None
  if cs is None or cs.unit_cell() is None or force_p1:
    cs, lo, margin = box_symmetry(hierarchy)
    offset = [margin - lo[i] for i in range(3)]
    print("no usable CRYST1: using a P1 box %s"
          % (" ".join("%.2f" % v for v in cs.unit_cell().parameters()[:3])))
  uc = cs.unit_cell()
  sg = cs.space_group()
  sgi = cs.space_group_info()

  # a centred description SHELX cannot express is not worth guessing at
  if sg.is_centric() and not sg.is_origin_centric():
    print("centric group with the origin off the centre: falling back to P1")
    cs = cs.customized_copy(space_group_info=sgtbx.space_group_info("P 1"))
    sg = cs.space_group()
    sgi = cs.space_group_info()

  # gather atoms first, so SFAC and UNIT can be written before them
  residues = []
  counts = {}
  n_atoms = 0
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          name = ag.resname.strip()
          if not keep_water and name in ("HOH", "DOD", "WAT"):
            continue
          atoms = []
          used = set()
          for a in ag.atoms():
            el = a.element.strip().capitalize()
            if not el:
              el = a.name.strip()[:1].upper()
            if el == "D":
              el = "H"
            site = a.xyz
            if offset is not None:
              site = tuple(site[i] + offset[i] for i in range(3))
            frac = uc.fractionalize(site)
            # PDB B, SHELX U
            u_iso = max(0.0, a.b/(8*math.pi*math.pi))
            atoms.append((sanitise_label(a.name, used), el, frac, u_iso,
                          a.occ))
            counts[el] = counts.get(el, 0) + 1
            n_atoms += 1
          if atoms:
            residues.append((name, chain.id.strip(), rg.resseq_as_int(),
                             rg.icode.strip(), atoms))
  if n_atoms == 0:
    raise RuntimeError("no atoms found in %s" % in_path)

  # SFAC order: C, H, then ascending Z
  elements = [e for e in PRIOR_ELEMENTS if e in counts]
  rest = sorted((tiny_pse.table(e).atomic_number(), e)
                for e in counts if e not in PRIOR_ELEMENTS)
  elements += [e for _, e in rest]
  sf_idx = dict((e, i + 1) for i, e in enumerate(elements))

  latt = 1 + "PIRFABC".find(sg.conventional_centring_type_symbol())
  if not sg.is_origin_centric():
    latt = -latt

  n_icode = 0
  with open(out_path, "w") as f:
    f.write("TITL %s in %s\n"
            % (os.path.basename(in_path), sgi.type().lookup_symbol()))
    f.write("REM converted by pdb2ins.py from %s\n" % in_path)
    f.write("CELL 0.71073 %.4f %.4f %.4f %.3f %.3f %.3f\n"
            % uc.parameters())
    f.write("ZERR %i 0. 0. 0. 0. 0. 0.\n" % sg.order_z())
    f.write("LATT %i\n" % latt)
    for i in range(sg.n_smx()):
      op = sg(0, 0, i)
      if op.is_unit_mx():
        continue
      f.write("SYMM %s\n" % op)
    f.write("SFAC %s\n" % " ".join(elements))
    f.write("UNIT %s\n" % " ".join("%d" % counts[e] for e in elements))
    f.write("\n")

    for (name, chain_id, seq, icode, atoms) in residues:
      if icode:
        n_icode += 1
      # RESI <class> [<chain>:]<number>, as TResidue::ToString writes it
      if chain_id:
        f.write("RESI %s %s:%d\n" % (name, chain_id, seq))
      else:
        f.write("RESI %s %d\n" % (name, seq))
      for (label, el, frac, u_iso, occ) in atoms:
        # occupancy fixed in the SHELX manner: 10 + q
        f.write("%-4s %2d %11.6f %11.6f %11.6f %11.5f %10.5f\n"
                % (label, sf_idx[el], frac[0], frac[1], frac[2],
                   10.0 + occ, u_iso))
    f.write("\nHKLF 4\nEND\n")

  print("%s -> %s" % (os.path.basename(in_path), out_path))
  print("  %d residues, %d atoms, elements %s"
        % (len(residues), n_atoms, " ".join(elements)))
  print("  cell %.2f %.2f %.2f  %s"
        % (uc.parameters()[0], uc.parameters()[1], uc.parameters()[2],
           sgi.type().lookup_symbol()))
  if n_icode:
    print("  WARNING: %d residues carry an insertion code, which SHELX cannot"
          " express. They keep their sequence number, so a chain may contain"
          " repeated numbers." % n_icode)
  return len(residues), n_atoms


def main(argv):
  if not argv:
    print(__doc__)
    return 1
  in_path = argv[0]
  out_path = None
  force_p1 = "--p1" in argv
  keep_water = "--no-water" not in argv
  for a in argv[1:]:
    if not a.startswith("--"):
      out_path = a
  if out_path is None:
    out_path = os.path.splitext(in_path)[0] + ".ins"
  convert(in_path, out_path, force_p1, keep_water)
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
