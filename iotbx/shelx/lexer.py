""" Parsing of ins/res files """
from __future__ import absolute_import, division, print_function

import re

from iotbx.shelx.errors import error as shelx_error
from iotbx.shelx import tokens

class command_stream(object):
  """ An ins/res file parsed as a stream of commands """

  commands_allowing_atom_names = { cmd: 1 for cmd in [
    'ANIS', 'BIND', 'BLOC', 'BOND', 'CHIV', 'CONF', 'CONN', 'DANG', 'DELU',
    'DFIX', 'EADP', 'EXYZ', 'FLAT', 'FREE', 'HFIX', 'HTAB', 'ISOR', 'MPLA',
    'NCSY', 'OMIT', 'RTAB', 'SADI', 'SAME', 'SIMU', 'RIGU'
  ]}

  shelx_commands = { cmd: 1 for cmd in [
    'ACTA', 'AFIX', 'BASF', 'BUMP', 'CELL', 'CGLS', 'DAMP', 'DEFS', 'DISP',
    'END' , 'EQIV', 'EXTI', 'FEND', 'FMAP', 'FRAG', 'FVAR', 'GRID', 'HKLF',
    'HOPE', 'L.S.', 'LATT', 'LAUE', 'LIST', 'MERG', 'MOLE', 'MORE', 'MOVE',
    'MUST', 'PART', 'PLAN', 'REM' , 'RESI', 'SFAC', 'SHEL', 'SIZE', 'SPEC',
    'STIR', 'SUMP', 'SWAT', 'SYMM', 'TEMP', 'TIME', 'TITL', 'TWIN', 'UNIT',
    'WGHT', 'WPDB', 'ZERR'
  ]}
  shelx_commands.update(commands_allowing_atom_names)


  def __init__(self, file=None, filename=None, close_when_deleted=True):
    assert [file, filename].count(None) == 1
    if file is None: file = open(filename)
    self.file = file
    self.close_when_deleted = close_when_deleted

  def __del__(self):
    if hasattr(self,'close_when_deleted') and self.close_when_deleted:
      if hasattr(self,'file') and hasattr(self.file,'close'):
        self.file.close()

  _cmd_pat = re.compile(r"""
    ^
    (?:
      (TITL | REM)
      \s?
      (.*)
      |
      (?:
        (?:
          ( %s )
          _
          (?:
            ([a-z] \S{0,3})
            |
            (\d{1,4})
            |
            \*
          )?
        )
        |
        ([a-z] \S{0,3})
      )
      ([^!=]*)
      (=?)
    )
    """ %(" | ".join(commands_allowing_atom_names.keys())),
        re.X | re.I)

  _continuation_pat = re.compile(r"""
    ^
    \s+
    ([^!=]*)
    (=?)
    """, re.X | re.I)

  _atom_name_pat = re.compile(r"""
    (?P<name>[a-z] [^_]{0,3})
    (?: _ \$ (?P<symmetry>\d+) )?
    (?: _ (?P<resnum>\d+))?
    (?: _ (?P<plusminus>[\+\-]))?
    |
    \$ (?P<element>[a-z]{1,2})
    |
    (?P<range> > | <)
    """, re.X | re.I)

  symm_space = re.compile(r"""
    (?<! , ) \s+
    """, re.X | re.I)

  _eqiv_pat = re.compile(r"""
    \$ (\d+) \s+ (.*) $
    """, re.X | re.I)

  _include_filename_pat = re.compile(r"""(\+)(.*)""")

  def __iter__(self):
    """
    Yields the commands in self.file as tuples:
      - either (cmd, args) where args is a tuple,
        e.g. ('CELL', (0.71073, 1, 2, 3, 90, 90, 90))
      - or (cmd, (residue_class, residue_number), args) for those commands
        which can be suffixed by a residue, e.g.
        ('HFIX', (parser.residue_number_tok, 1), (23, ))
    Notes:
      - For atoms, cmd is '__ATOM__' and the name is the first argument
      - For so-called Q-peaks, cmd is '__Q_PEAK__' and the number is the
        first argument
      - In args, floating point items are reported as is whereas any
        string comes as (type, value) where type is one of the class
        constants defined just above this method (actually, for
        the residue_number_tok, the value is an integer stored as a float).
      - The atomic refinable variables are not decoded: this is the job of
        a parser, not of a lexer.
    """
    continued = False
    for i, li in enumerate(self.file):
      if not li: continue
      if continued:
        m = self._continuation_pat.search(li)
        if m is None:
          raise shelx_error("illegal continuation line error", i)
        cont_args, continued = m.groups()
        arguments.extend(cont_args.split())
      else:
        m = self._cmd_pat.search(li)
        if m is None:
          if li[0].isspace(): continue
          else:
            m = self._include_filename_pat.search(li)
            if m is not None:
              cmd, filename = m.group(1, 2)
              yield (cmd, filename), i
              continue
          raise shelx_error("illegal command or atom name error", i)
        cmd = m.group(1) or m.group(6)
        if cmd:
          cmd = cmd.upper()
          cmd_residue = None
        else:
          (cmd, cmd_residue_class, cmd_residue_number) = m.group(3,4,5)
          cmd = cmd.upper()
          if cmd_residue_class:
            cmd_residue = (tokens.residue_class_tok, cmd_residue_class.upper())
          elif cmd_residue_number:
            cmd_residue = (tokens.residue_number_tok, int(cmd_residue_number))
          else:
            cmd_residue = (tokens.all_residues_tok,)
        args = m.group(2) or m.group(7) or ""
        continued = m.group(8)
        arguments = args.split()
      if not continued:
        result = self._parse_special_cases(cmd, args, arguments, i, li)
        if result is None:
          result = self._parse_general_case(cmd, cmd_residue, arguments, i, li)
        yield result, i
        if cmd == 'HKLF': break

  def _parse_special_cases(self, cmd, args, arguments, i, li):
    if cmd in ('TITL', 'REM'):
      args = args.strip()
      if args: arg_tuple = (args,)
      else: arg_tuple = ()
      return (cmd, arg_tuple)
    if cmd == 'SYMM':
      if args is None:
        raise shelx_error("illegal argument '%s'", i, args)
      return (cmd, (self.symm_space.sub('', args).upper(),))
    if cmd == 'EQIV':
      m = self._eqiv_pat.search(args)
      if m is None:
        raise shelx_error("illegal argument '%s'", i, args)
      try:
        idx = int(m.group(1))
      except ValueError:
        raise shelx_error("illegal argument '%s'", i, m.group(1))
      return (cmd, (idx, self.symm_space.sub('', m.group(2)).upper()))
    if cmd == 'SFAC':
      sfac_args = []
      for e in arguments:
        try:
          sfac_args.append(float(e))
        except ValueError:
          sfac_args.append(e)
      return (cmd, tuple(sfac_args))
    return None

  def _parse_general_case(self, cmd, cmd_residue, arguments, i, li):
    toks = []
    for j,arg in enumerate(arguments):
      try:
        toks.append(float(arg))
      except ValueError:
        if cmd == 'RESI':
          toks.append(arg)
        else:
          m = self._atom_name_pat.match(arg)
          if m is None:
            raise shelx_error("illegal argument '%s'", i, arg)
          name, symmetry, resnum, plus_minus = m.group(
            'name', 'symmetry', 'resnum', 'plusminus')
          if name is not None:
            if resnum is not None: resnum = int(resnum)
            if symmetry is not None: symmetry = int(symmetry)
            toks.append(tokens.atomname_token(name=name,
                                              symmetry=symmetry,
                                              residue_number=resnum,
                                              plus_minus=plus_minus))
            continue
          element, resnum = m.group('element', 'resnum')
          if element is not None:
            toks.append(tokens.element_token(element=element,
                                             residue_number=resnum))
            continue

          if m.group('range') == '>':
            toks.append((tokens.forward_range_tok,))
          elif m.group('range') == '<':
            toks.append((tokens.backward_range_tok,))
    if cmd_residue:
      return (cmd, cmd_residue, tuple(toks))
    else:
      if cmd not in self.shelx_commands:
        if cmd[0] == 'Q':
          toks = (int(cmd[1:]),) + tuple(toks)
          cmd = '__Q_PEAK__'
        else:
          toks = (cmd,) + tuple(toks)
          cmd = '__ATOM__'
      else:
        toks = tuple(toks)
      return (cmd, toks)
