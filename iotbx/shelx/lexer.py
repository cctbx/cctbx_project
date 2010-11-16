""" Parsing of ins/res files """

import re

from iotbx.shelx.errors import error as shelx_error

class command_stream(object):
  """ An ins/res file parsed as a stream of commands """

  shelx_commands = dict([ (cmd, 1) for cmd in [
    'ACTA', 'AFIX', 'ANIS', 'BASF', 'BIND', 'BLOC', 'BOND', 'BUMP', 'CELL',
    'CGLS', 'CHIV', 'CONF', 'CONN', 'DAMP', 'DANG', 'DEFS', 'DELU', 'DFIX',
    'DISP', 'EADP', 'END ', 'EQIV', 'EXTI', 'EXYZ', 'FEND', 'FLAT', 'FMAP',
    'FRAG', 'FREE', 'FVAR', 'GRID', 'HFIX', 'HKLF', 'HOPE', 'HTAB', 'ISOR',
    'L.S.', 'LATT', 'LAUE', 'LIST', 'MERG', 'MOLE', 'MORE', 'MOVE', 'MPLA',
    'MUST', 'NCSY', 'OMIT', 'PART', 'PLAN', 'REM ', 'RESI', 'RTAB', 'SADI',
    'SAME', 'SFAC', 'SHEL', 'SIMU', 'SIZE', 'SPEC', 'STIR', 'SUMP', 'SWAT',
    'SYMM', 'TEMP', 'TIME', 'TITL', 'TWIN', 'UNIT', 'WGHT', 'WPDB', 'ZERR'
  ]])

  def __init__(self, file=None, filename=None):
    assert [file, filename].count(None) == 1
    if file is None: file = open(filename)
    self.file = file

  _cmd_pat = re.compile("""
    ^
    (?:
      (?:
        ( HFIX | EXYZ | EADP | OMIT | CONN |
          SAME | CHIV | DELU | SIMU | ISOR |
          BLOC | BOND | CONF | MPLA | RTAB )
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
    """, re.X | re.I)

  _continuation_pat = re.compile("""
    ^
    \s+
    ([^!=]*)
    (=?)
    """, re.X | re.I)

  _atom_name_pat = re.compile("""
    ([a-z] [^_]{0,3})
    (?: _ \$ (\d+) )?
    |
    \$ ([a-z]{1,2})
    |
    (> | <)
    """, re.X | re.I)

  symm_space = re.compile("""
    (?<! , ) \s+
    """, re.X | re.I)

  _eqiv_pat = re.compile("""
    \$ (\d+) \s+ (.*) $
    """, re.X | re.I)

  (atom_tok, element_tok,
   forward_range_tok, backward_range_tok,
   residue_class_tok, residue_number_tok, all_residues_tok) = xrange(7)

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
          raise shelx_error("illegal command or atom name error", i)
        cmd = m.group(4)
        if cmd:
          cmd = cmd.upper()
          cmd_residue = None
        else:
          (cmd, cmd_residue_class, cmd_residue_number) = m.group(1,2,3)
          cmd = cmd.upper()
          if cmd_residue_class:
            cmd_residue = (self.residue_class_tok, cmd_residue_class.upper())
          elif cmd_residue_number:
            cmd_residue = (self.residue_number_tok, int(cmd_residue_number))
          else:
            cmd_residue = (self.all_residues_tok,)
        args, continued = m.group(5,6)
        arguments = args.split()
      if not continued:
        result = self._parse_special_cases(cmd, args, i, li)
        if result is None:
          result = self._parse_general_case(cmd, cmd_residue, arguments, i, li)
        yield result, i
        if cmd == 'HKLF': break

  def _parse_special_cases(self, cmd, args, i, li):
    if cmd in ('TITL', 'REM'):
      return (cmd, (args.strip(),))
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
      return (cmd, tuple([ e.upper() for e in args.split() ]))
    return None

  def _parse_general_case(self, cmd, cmd_residue, arguments, i, li):
    tokens = []
    for j,arg in enumerate(arguments):
      try:
        tokens.append(float(arg))
      except ValueError:
        if cmd == 'RESI':
          tokens.append(arg)
        else:
          m = self._atom_name_pat.match(arg)
          if m is None:
            raise shelx_error("illegal argument '%s'", i, arg)
          name, symmetry = m.group(1,2)
          if name is not None:
            if symmetry is not None: symmetry = int(symmetry)
            tokens.append((self.atom_tok, name, symmetry))
            continue
          element = m.group(3)
          if element is not None:
            tokens.append((self.element_tok, element))
            continue
          if arg == '>':
            tokens.append((self.forward_range_tok,))
          else:
            tokens.append((self.backward_range_tok,))
    if cmd_residue:
      return (cmd, cmd_residue, tuple(tokens))
    else:
      if len(cmd) < 4 or cmd not in self.shelx_commands:
        tokens = (cmd,) + tuple(tokens)
        cmd = '__ATOM__'
      else:
        tokens = tuple(tokens)
      return (cmd, tokens)
