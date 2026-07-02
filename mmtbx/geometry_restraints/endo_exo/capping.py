"""Hydrogen cap placement at covalent boundary sites."""

from __future__ import absolute_import, division, print_function

from cctbx.array_family import flex


class HydrogenCapper:
  """Place hydrogen cap atoms at covalent boundary sites.

  Parameters
  ----------
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, log=None):
    self.log = log

  def cap_atom(self, anchor, cap):
    """Move *cap* to a hydrogen position 1.1 A along the anchor->cap vector.

    Parameters
    ----------
    anchor : iotbx.pdb.hierarchy.atom or None
    cap : iotbx.pdb.hierarchy.atom or None
        Both may be ``None`` (no-op).
    """
    if anchor is None or cap is None:
      return
    v = flex.double(cap.xyz) - flex.double(anchor.xyz)
    v_norm = v.norm()
    assert v_norm > 1e-6, "anchor and cap must be distinct atoms"
    u = v / v_norm
    cap.set_element('H')
    # cap.set_name(('H' + cap.name.strip()).rjust(4))
    cap.set_xyz(tuple(flex.double(anchor.xyz) + u * 1.1))
    print('Capped atom:', file=self.log)
    print('  ' + anchor.format_atom_record().rstrip(), file=self.log)
    print('  ' + cap.format_atom_record().rstrip(), file=self.log)
