"""Amino-acid sidechain and free-terminus net-charge estimation."""

from __future__ import absolute_import, division, print_function

from collections import defaultdict


STANDARD_RESIDUE_CHARGES = {
  'ALA':  0,
  'ARG':  1,
  'ASN':  0,
  'ASP': -1,
  'CYS':  0,
  'GLN':  0,
  'GLU': -1,
  'GLY':  0,
  'HIS':  None,  # can be +1, 0 or -1 depending on protonation state
  'ILE':  0,
  'LEU':  0,
  'LYS':  1,
  'MET':  0,
  'PHE':  0,
  'PRO':  0,
  'SER':  0,
  'THR':  0,
  'TRP':  0,
  'TYR':  0,
  'VAL':  0,
}

# formal_charge is the residue's charge in its standard RCSB CCD reference
# form (ASP/GLU/CYS/TYR neutral; LYS/ARG/HIS+ protonated).  protonation_sites
# lists every canonical H name that, when present, contributes +1 to the
# reference charge -- one canonical name per chemically distinct H slot.
# Net sidechain charge = formal_charge - len(canonical_sites) + n_present.
CHARGED_SIDECHAINS = {
  'ASP': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OD1', 'OD2'},
    'protonation_sites': {'HD2', 'DD2'},
  },
  'CYS': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'SG'},
    'protonation_sites': {'HG', 'DG'},
  },
  'GLU': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OE1', 'OE2'},
    'protonation_sites': {'HE2', 'DE2'},
  },
  'LYS': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'NZ'},
    'protonation_sites': {'HZ1', 'HZ2', 'HZ3', 'DZ1', 'DZ2', 'DZ3'},
  },
  'ARG': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'NE', 'NH1', 'NH2', 'CZ'},
    'protonation_sites': {'HH11', 'HH12', 'HH21', 'HH22', 'DH11', 'DH12', 'DH21', 'DH22'},
  },
  'HIS': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'ND1', 'NE2'},
    'protonation_sites': {'HD1', 'HE2', 'DD1', 'DE2'},
  },
  'TYR': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OH'},
    'protonation_sites': {'HH', 'DH'},
  },
}


class ChargeEstimator:
  """Estimate the net amino-acid sidechain charge of a model.

  Walks the model's residues, looks up each standard amino acid's expected
  charge from ``CHARGED_SIDECHAINS`` / ``STANDARD_RESIDUE_CHARGES`` based on
  which protonation-site H atoms are actually present, and (optionally)
  adds a contribution for free peptide termini.  The model may be a full
  structure, a truncated QM region, or any subset thereof.

  Parameters
  ----------
  include_terminal_charges : bool, optional
      If ``True``, detect free peptide termini and add their charges.
      Default is ``False``.
  n_terminus_charge : int, optional
      Charge assigned to each detected free N-terminus.  Default is ``1``.
  c_terminus_charge : int, optional
      Charge assigned to each detected free C-terminus.  Default is ``-1``.
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, include_terminal_charges=False,
               n_terminus_charge=1, c_terminus_charge=-1, log=None):
    self.include_terminal_charges = include_terminal_charges
    self.n_terminus_charge = n_terminus_charge
    self.c_terminus_charge = c_terminus_charge
    self.log = log

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  def calculate(self, model):
    """Estimate the net amino-acid sidechain charge of *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    dict
        Keys: ``total_charge``, ``sidechain_charge``, ``terminal_charge``,
        ``contributors``, ``residue_counts``, ``unknown_residues``,
        ``residue_count_total``.
    """
    residue_data = self._collect_residue_data(model.get_hierarchy().atoms())

    total_sidechain_charge = 0
    total_terminal_charge = 0
    residue_counts = defaultdict(int)
    unknown_residues = set()
    residue_contributions = {}

    for key, data in residue_data.items():
      resname = data['resname']
      if resname not in STANDARD_RESIDUE_CHARGES:
        unknown_residues.add(resname)
        continue
      residue_counts[resname] += 1
      sidechain_charge = self._sidechain_charge(resname, data['atom_names'])
      total_sidechain_charge += sidechain_charge
      residue_contributions[key] = {
        'resname': resname,
        'sidechain_charge': sidechain_charge,
        'terminal_charge': 0,
      }

    if self.include_terminal_charges:
      total_terminal_charge = self._add_terminal_charges(
        residue_data, residue_contributions, model
      )

    contributors = self._build_contributors_list(residue_contributions)
    total_charge = total_sidechain_charge + total_terminal_charge

    return {
      'total_charge': total_charge,
      'sidechain_charge': total_sidechain_charge,
      'terminal_charge': total_terminal_charge,
      'contributors': contributors,
      'residue_counts': dict(sorted(residue_counts.items())),
      'unknown_residues': sorted(unknown_residues),
      'residue_count_total': len(residue_data),
    }

  def show(self, seed_index, charge_summary, out=None):
    """Write a human-readable charge summary to *out*.

    No-op when *out* is ``None``.

    Parameters
    ----------
    seed_index : int
    charge_summary : dict
        As returned by :meth:`calculate`.
    out : file-like or None, optional
    """
    if out is None:
      return
    fmt = self._fmt_signed
    print(
      f'Estimated amino-acid net charge (seed {seed_index}): '
      f'{fmt(charge_summary["total_charge"])}',
      file=out,
    )
    print(f'  sidechain contribution: {fmt(charge_summary["sidechain_charge"])}',
          file=out)
    if self.include_terminal_charges:
      print(f'  terminal contribution : {fmt(charge_summary["terminal_charge"])}',
            file=out)
    print(f'  residues counted      : {charge_summary["residue_count_total"]}',
          file=out)
    if charge_summary['residue_counts']:
      counts_str = ', '.join(
        f'{resname}:{count}'
        for resname, count in charge_summary['residue_counts'].items()
      )
      print(f'  residue composition   : {counts_str}', file=out)
    if charge_summary['unknown_residues']:
      print(
        '  skipped non-standard residues: ' +
        ', '.join(charge_summary['unknown_residues']),
        file=out,
      )
    if charge_summary['contributors']:
      print('  charge-contributing residues:', file=out)
      for entry in charge_summary['contributors']:
        residue_id = f'chain {entry["chain_id"]} resseq {entry["resseq"]}'
        if entry['icode']:
          residue_id += f' icode {entry["icode"]}'
        if entry['altloc']:
          residue_id += f' altloc {entry["altloc"]}'
        print(
          f'    {entry["resname"]} {residue_id}: '
          f'sidechain: {fmt(entry["sidechain_charge"])}, '
          f'terminal: {fmt(entry["terminal_charge"])}; '
          f'total: {fmt(entry["total_charge"])}',
          file=out,
        )

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  @staticmethod
  def _fmt_signed(x):
    return f'{x:+}' if x else '0'

  @staticmethod
  def _collect_residue_data(atoms):
    """Group atoms by residue and collect atom names and backbone i_seqs.

    Parameters
    ----------
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    dict
        Keyed by ``(chain_id, resseq, icode, altloc)`` tuples.
    """
    residue_data = {}
    for atom in atoms:
      i_seq = atom.i_seq
      atom_group = atom.parent()
      residue_group = atom_group.parent()
      chain = residue_group.parent()
      resname = atom_group.resname.strip().upper()

      key = (
        chain.id.strip(),
        residue_group.resseq.strip(),
        residue_group.icode.strip(),
        atom_group.altloc.strip(),
      )
      if key not in residue_data:
        residue_data[key] = {
          'resname': resname,
          'atom_names': set(),
          'n_i_seqs': set(),
          'c_i_seqs': set(),
        }

      atom_name = atom.name.strip().upper()
      residue_data[key]['atom_names'].add(atom_name)
      if atom_name == 'N':
        residue_data[key]['n_i_seqs'].add(i_seq)
      elif atom_name == 'C':
        residue_data[key]['c_i_seqs'].add(i_seq)

    return residue_data

  @staticmethod
  def _is_charged_sidechain_present(resname, atom_names):
    """Return ``True`` if the charge-bearing heavy atom(s) of *resname* are
    in *atom_names*.

    Residues without a ``charged_heavy_atoms`` entry in
    ``CHARGED_SIDECHAINS`` are considered always present (return ``True``);
    otherwise the residue's sidechain is treated as truncated/absent when
    none of its charged heavy atoms appear.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case).
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    bool
    """
    charged_atoms = CHARGED_SIDECHAINS.get(resname, {}).get(
      'charged_heavy_atoms', set()
    )
    if not charged_atoms:
      return True
    return bool(charged_atoms.intersection(atom_names))

  @staticmethod
  def _calculate_side_chain_charge(resname, atom_names):
    """Compute the sidechain net charge from the present protonation-site H atoms.

    Starts from ``formal_charge`` (the residue's charge in its CCD
    reference protonation state) and adjusts for missing or extra
    protonation-site hydrogens.  Deuterium names (``D...``) are folded
    onto their canonical hydrogen names so D/H labelling does not affect
    the count.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case); must be a key in
        ``CHARGED_SIDECHAINS``.
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    int
        Net sidechain charge.
    """
    entry = CHARGED_SIDECHAINS.get(resname, {})
    formal_charge = entry.get('formal_charge', 0)
    protonation_sites = entry.get('protonation_sites', set())

    def _canonical_h_site(name):
      return 'H' + name[1:] if name.startswith('D') else name

    sites = {_canonical_h_site(n) for n in protonation_sites}
    present = {_canonical_h_site(n) for n in atom_names}
    return formal_charge - len(sites) + len(sites.intersection(present))

  def _sidechain_charge(self, resname, atom_names):
    """Return the net sidechain charge for *resname* given its present atoms.

    Returns ``0`` if the residue's charge-bearing heavy atoms are absent
    (truncated sidechain).  Otherwise, residues listed in
    ``CHARGED_SIDECHAINS`` get a protonation-aware charge; remaining
    residues fall back to ``STANDARD_RESIDUE_CHARGES``.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case).
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    int
    """
    if not self._is_charged_sidechain_present(resname, atom_names):
      return 0
    if resname in CHARGED_SIDECHAINS:
      return self._calculate_side_chain_charge(resname, atom_names)
    return STANDARD_RESIDUE_CHARGES.get(resname, 0)

  def _add_terminal_charges(self, residue_data, residue_contributions, model):
    """Detect free termini and accumulate terminal charges.

    A residue is treated as having a free N-terminus if it is the first
    residue group in its chain (within *model*), and likewise a free
    C-terminus if it is the last.

    Parameters
    ----------
    residue_data : dict
    residue_contributions : dict
        Modified in-place.
    model : mmtbx.model.manager

    Returns
    -------
    int
        Total terminal charge contribution.
    """
    first_rgs, last_rgs = self._chain_terminus_keys(model)
    total = 0
    for key, data in residue_data.items():
      if key not in residue_contributions:
        continue
      rg_key = key[:3]  # (chain_id, resseq, icode); altloc is dropped

      if data['n_i_seqs'] and rg_key in first_rgs:
        residue_contributions[key]['terminal_charge'] += self.n_terminus_charge
        total += self.n_terminus_charge
      if data['c_i_seqs'] and rg_key in last_rgs:
        residue_contributions[key]['terminal_charge'] += self.c_terminus_charge
        total += self.c_terminus_charge

    return total

  @staticmethod
  def _chain_terminus_keys(model):
    """Return ``(first, last)`` sets of residue-group keys for chain termini.

    Each key is ``(chain_id, resseq, icode)`` for the first / last residue
    group in its chain within *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    first : set of tuple
    last : set of tuple
    """
    first = set()
    last = set()
    for chain in model.get_hierarchy().chains():
      rgs = list(chain.residue_groups())
      if not rgs:
        continue
      chain_id = chain.id.strip()
      first.add((chain_id, rgs[0].resseq.strip(), rgs[0].icode.strip()))
      last.add((chain_id, rgs[-1].resseq.strip(), rgs[-1].icode.strip()))
    return first, last

  @staticmethod
  def _build_contributors_list(residue_contributions):
    """Flatten *residue_contributions* into a sorted list of nonzero entries.

    Residues whose total charge (sidechain + terminal) is zero are
    omitted from the output.

    Parameters
    ----------
    residue_contributions : dict
        Keyed by ``(chain_id, resseq, icode, altloc)``; values carry
        ``resname``, ``sidechain_charge``, and ``terminal_charge``.

    Returns
    -------
    list of dict
        One entry per nonzero-charge residue, sorted by key.  Each entry
        has keys ``chain_id``, ``resseq``, ``icode``, ``altloc``,
        ``resname``, ``sidechain_charge``, ``terminal_charge``,
        ``total_charge``.
    """
    contributors = []
    for key, contribution in sorted(residue_contributions.items()):
      chain_id, resseq, icode, altloc = key
      sc = contribution['sidechain_charge']
      tc = contribution['terminal_charge']
      total = sc + tc
      if abs(total) == 0:
        continue
      contributors.append({
        'chain_id': chain_id,
        'resseq': resseq,
        'icode': icode,
        'altloc': altloc,
        'resname': contribution['resname'],
        'sidechain_charge': sc,
        'terminal_charge': tc,
        'total_charge': total,
      })
    return contributors
