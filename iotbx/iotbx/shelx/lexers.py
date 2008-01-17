""" Lexing of ins/res files """

from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import adptbx

import scitbx.math

from libtbx import forward_compatibility
from libtbx import adopt_init_args

class crystal_symmetry_lexer(object):

  def __init__(self, parser):
    self.parser = parser

  def filtered_commands(self):
    unit_cell = None
    space_group = sgtbx.space_group()
    for command in self.parser.commands():
      cmd, args = command[0], command[-1]
      if cmd == 'CELL':
        assert unit_cell is None
        unit_cell = uctbx.unit_cell(args[1:])
      elif cmd == 'LATT':
        assert len(args) == 1
        n = int(args[0])
        if n > 0:
          space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
        z = "*PIRFABC"[abs(n)]
        space_group.expand_conventional_centring_type(z)
      elif cmd == 'SYMM':
        assert len(args) == 1
        s = sgtbx.rt_mx(args[0])
        space_group.expand_smx(s)
      else:
        if cmd == 'SFAC':
          assert unit_cell is not None
          self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                                   space_group=space_group)
        yield command

  def lex(self):
    for command in self.filtered_commands():
      if command[0] == 'SFAC': break


class atom_lexer(crystal_symmetry_lexer):

  def __init__(self,
               parser,
               set_grad_flags=False,
               min_distance_sym_equiv=0.5):
    adopt_init_args(self, locals())

  def filtered_commands(self):
    self.structure = None
    for command in super(atom_lexer, self).filtered_commands():
      cmd, args = command[0], command[-1]
      if cmd == 'SFAC':
        self.structure = xray.structure(
          special_position_settings=crystal.special_position_settings(
            crystal_symmetry=self.crystal_symmetry,
            min_distance_sym_equiv=self.min_distance_sym_equiv))
        self.label_for_sfac = ('*',) + args # (a) working around
                                            #     ShelXL 1-based indexing
      elif cmd == 'FVAR':
        self.overall_scale = args[0]
        self.free_variable = args # (b) actually args[1:] but working around
                                  #     ShelXL 1-based indexing
      elif cmd == 'PART' and len(args) == 2:
        raise NotImplementedError
      else:
        atom = False
        if len(cmd) < 4:
          atom = True
        elif len(args) in (6, 11):
          atom = True
        if not atom:
          yield command
          continue
        try:
          self.lex_atom(cmd, args)
        except atom_lexer.not_an_atom:
          yield command
          continue

  class not_an_atom(RuntimeError): pass

  def lex_atom(self, cmd, args):
    n = int(args[0])
    n_vars = len(args) - 1
    if n_vars == 5:
      values, behaviours = self.decode_variables(args[1:], u_iso_idx=n_vars-1)
      u = values[-1]
      self.previous_boundable_u_eq = u
      isotropic = True
    else:
      unit_cell = self.crystal_symmetry.unit_cell()
      values, behaviours = self.decode_variables(
        args[1:-3] + (args[-1], args[-2], args[-3]),
        u_iso_idx=None)
      u = adptbx.u_cif_as_u_star(unit_cell, values[-6:])
      self.previous_boundable_u_eq = adptbx.u_star_as_u_iso(unit_cell, u)
      isotropic = False
    site = values[0:3]
    occ = values[3]
    scattering_type = eltbx.xray_scattering.get_standard_label(
      self.label_for_sfac[n], # works thank to (a)
      exact=True)
    scatterer = xray.scatterer(
      label           = cmd,
      site            = site,
      occupancy       = occ,
      u               = u,
      scattering_type = scattering_type)
    if self.set_grad_flags:
      f = scatterer.flags
      if atom_lexer.fixed not in behaviours[0:3]:
        f.set_grad_site(True)
      if behaviours[3] != atom_lexer.fixed:
        f.set_grad_occupancy(True)
      if isotropic:
        if behaviours[4] != atom_lexer.fixed:
          f.set_grad_u_iso(True)
      else:
        if atom_lexer.fixed not in behaviours[-6:]:
          f.set_grad_u_aniso(True)
    self.structure.add_scatterer(scatterer)

  (freely_refined, fixed, bound_to_previous_u_eq,
   bound_to_free_variable) = xrange(-1,3)

  def decode_variables(self, coded_variables, u_iso_idx=None):
    # this method works thanks to (b) above
    values = []
    behaviours = []
    for i,coded_variable in enumerate(coded_variables):
      try:
        m,p = scitbx.math.divmod(coded_variable, 10)
      except ArgumentError:
        raise atom_lexer.not_an_atom
      if m <= -2:
        # p*(fv_{-m} - 1)
        m = -m-1
        values.append( p*(self.free_variable[m] - 1) )
        behaviours.append(m)
      elif m == 0:
        if i == u_iso_idx and p < -0.5:
          # p * (U_eq of the previous atom not constrained in this way)
          values.append( -p*self.previous_boundable_u_eq )
          behaviours.append(atom_lexer.bound_to_previous_u_eq)
        else:
          # p (free to refine)
          values.append(p)
          behaviours.append(atom_lexer.freely_refined)
      elif m == 1:
        # p (fixed variable)
        values.append(p)
        behaviours.append(atom_lexer.fixed)
      elif m >= 2:
        # p*fv_m
        m = m-1
        values.append(p*self.free_variable[m])
        behaviours.append(m)
      else:
        # m == -1
        # undocumented, rather pathological case
        # but I carefully checked that ShelXL does indeed behave so!
        values.append(0)
        behaviours.append(atom_lexer.fixed)
    return values, behaviours

  def lex(self):
    for command in self.filtered_commands(): pass
