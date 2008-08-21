""" Lexing of ins/res files """

from __future__ import generators

from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import adptbx

import scitbx.math

from libtbx import forward_compatibility
from libtbx import adopt_init_args
import libtbx.load_env

from iotbx.shelx import util
from iotbx.shelx.errors import error as shelx_error

if (libtbx.env.dist_path("smtbx", default=None) is None):
  smtbx = None
else:
  import smtbx
  import smtbx.refinement.constraints as smtbx_constraints


class parser(object):

  def __init__(self, command_stream, builder=None):
    self.command_stream = command_stream
    self.builder = builder

  def parse(self):
    for command, line in self.filtered_commands(): pass


class crystal_symmetry_parser(parser):
  """ A parser pulling out the crystal symmetry info from a command stream """

  def filtered_commands(self):
    """ Yields those command in self.command_stream
        that this parser is not concerned with. On the contrary,
        CELL, LATT, SYMM are swallowed.
        The resulting info is available in self.crystal_symmetry
    """
    unit_cell = None
    space_group = sgtbx.space_group()
    for command, line in self.command_stream:
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
          self.builder.make_crystal_symmetry(unit_cell=unit_cell,
                                             space_group=space_group)
        yield command, line

  def parse(self):
    for command, line in self.filtered_commands():
      if command[0] == 'SFAC': break


class variable_decoder(util.behaviour_of_variable):

  def decode_variables(self, coded_variables, u_iso_idx=None):
    values = []
    behaviours = []
    for i,coded_variable in enumerate(coded_variables):
      try:
        m,p = scitbx.math.divmod(coded_variable, 10)
      except ArgumentError:
        raise shelx_error("%i-th scatterer parameter '%s'",
                          self.line, i, coded_variable)
      if m <= -2:
        # p*(fv_{-m} - 1)
        m = -m-1 # indexing thanks to (b) above
        values.append( p*(self.free_variable[m] - 1) )
        behaviours.append((self.p_times_fv_minus_1, p, m))
      elif m == 0:
        if i == u_iso_idx and p < -0.5:
          # p * (U_eq of the previous atom not constrained in this way)
          scatt, scatt_idx = self.builder.scatterer_to_bind_u_eq_to
          u_iso = scatt.u_eq(self.builder.crystal_symmetry.unit_cell())
          values.append( -p*u_iso )
          behaviours.append((self.p_times_previous_u_eq, scatt_idx))
        else:
          # p (free to refine)
          values.append(p)
          behaviours.append(self.free)
      elif m == 1:
        # p (fixed variable)
        values.append(p)
        behaviours.append(self.fixed)
      elif m >= 2:
        # p*fv_m
        m = m-1 # indexing thanks to (b) above
        values.append(p*self.free_variable[m])
        behaviours.append((self.p_times_fv, p, m))
      else:
        # m == -1
        # undocumented, rather pathological case
        # but I carefully checked that ShelXL does indeed behave so!
        values.append(0)
        behaviours.append(self.fixed)
    return values, behaviours


class atom_parser(parser, variable_decoder):
  """ A parser pulling out the scatterer info from a command stream """

  def filtered_commands(self):
    self.label_for_sfac = None
    scatterer_index = 0
    for command, line in self.command_stream:
      self.line = line
      cmd, args = command[0], command[-1]
      if cmd == 'SFAC':
        self.builder.make_structure()
        self.label_for_sfac = ('*',) + args # (a) working around
                                            #     ShelXL 1-based indexing
      elif cmd == 'FVAR':
        self.overall_scale = args[0]
        self.free_variable = args # (b) ShelXL indexes into the whole array
      elif cmd == 'PART' and len(args) == 2:
        raise NotImplementedError
      elif cmd == '__ATOM__':
        if self.label_for_sfac is None:
          raise shelx_error("missing sfac", self.line)
        scatterer, behaviour_of_variable = self.lex_scatterer(
          args, scatterer_index)
        self.builder.add_scatterer(scatterer, behaviour_of_variable)
        scatterer_index += 1
      else:
        yield command, line

  def lex_scatterer(self, args, scatterer_index):
    name = args[0]
    n = int(args[1])
    n_vars = len(args) - 2
    if n_vars == 5:
      values, behaviours = self.decode_variables(
        args[2:],
        u_iso_idx=n_vars-1)
      u = values[-1]
      isotropic = True
    elif n_vars == 10:
      unit_cell = self.builder.crystal_symmetry.unit_cell()
      values, behaviours = self.decode_variables(
        args[2:-3] + (args[-1], args[-2], args[-3]),
        u_iso_idx=None)
      u = adptbx.u_cif_as_u_star(unit_cell, values[-6:])
      isotropic = False
    else:
      raise shelx_error("wrong number of parameters for scatterer",
                        self.line)
    site = values[0:3]
    occ = values[3]
    scattering_type = eltbx.xray_scattering.get_standard_label(
      self.label_for_sfac[n], # works thank to (a)
      exact=True)
    scatterer = xray.scatterer(
      label           = name,
      site            = site,
      occupancy       = occ,
      u               = u,
      scattering_type = scattering_type)
    if not isotropic or behaviours[-1] != self.p_times_previous_u_eq:
      self.builder.scatterer_to_bind_u_eq_to = (scatterer, scatterer_index)
    return scatterer, behaviours


if smtbx is not None:

  class afix_parser(parser):
    """ It must be before an atom parser """

    constraint = {
      1:  (smtbx_constraints.tertiary_CH, 1),
      2:  (smtbx_constraints.secondary_CH2, 2),
      3:  (smtbx_constraints.staggered_terminal_tetrahedral_XHn, 3),
      4:  (smtbx_constraints.aromatic_CH_or_amide_NH, 1),
      8:  (smtbx_constraints.staggered_terminal_tetrahedral_XHn, 1),
      9:  (smtbx_constraints.terminal_trihedral_XH2, 2),
      13: (smtbx_constraints.terminal_tetrahedral_XHn, 3),
      14: (smtbx_constraints.terminal_tetrahedral_XHn, 1),
      15: (smtbx_constraints.polyhedral_BH, 1),
      16: (smtbx_constraints.acetylenic_CH, 1),
      }

    def filtered_commands(self):
      active_afix = False
      for command, line in self.command_stream:
        cmd, args = command[0], command[-1]
        if cmd in ('AFIX', 'HKLF'):
          if cmd == 'AFIX' and not args:
            raise shelx_error("too few arguments", line)
          if active_afix:
            if n_afixed != n_expected_afixed:
              raise shelx_error("wrong number of afixed atoms", line)
            self.builder.end_afix()
            active_afix = False
            n_afixed = 0
          if cmd == 'HKLF':
            yield command, line
            continue
          mn = args[0]
          if mn == 0: continue
          m,n = divmod(mn, 10)
          d, sof, u = (None,)*3
          params = args[1:]
          if not params: pass
          elif len(params) == 1: d = params[0]
          elif len(params) == 3: d, sof = params[0:]
          elif len(params) == 4: d, sof, u = params[0:]
          else: raise shelx_error("too many arguments", line)
          info = self.constraint.get(m)
          if info is not None:
            constraint_type, n_expected_afixed = info
            kwds = {}
            if d is not None: kwds['bond_length'] = d
            if (n in (7,8)
                and not constraint_type.__name__.startswith('staggered')):
              kwds['rotating'] = True
            self.builder.start_afix(constraint_type, kwds)
            active_afix = True
            n_afixed = 0
        elif cmd == '__ATOM__':
          if active_afix:
            if sof is not None: args[4] = sof
            if u is not None and len(args) == 6: args[-1] = u
            n_afixed += 1
          yield command, line
        else:
          yield command, line
      self.builder.finish()
