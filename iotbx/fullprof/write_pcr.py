from __future__ import absolute_import, division, print_function
from cctbx.eltbx import wavelengths
from six.moves import range

def _make_phase_block(phase, number=1, name="", scale_down=1.0):
  """Create a pcr phase block skelleton with placeholder strings for different
  refinement options.

  (for internal use)

  :param phase: a crystal structucture
  :type phase: cctbx.xray_structure
  :param number: the number of the phase
  :type number: integer
  :param name: a title to be used for the phase in the pcr file
  :type name: string
  :param scale_down: factor to divide intensities by (to avoid overflows)
  :type scale_down: float

  :returns: the pcr phase block skelleton as a string
  :rtype: string
  """
  phase.make_scatterer_labels_shelx_compatible_in_place()
  if name =="":
    name = "Phase_{0}".format(number)
  scatt = phase.scatterers()
  symm = phase.crystal_symmetry()
  tmp = """\
!-------------------------------------------------------------------------------
!  Data for PHASE number:  {number}  ==> Current R_Bragg for Pattern#  {number}:
!-------------------------------------------------------------------------------
 {name:s}
!
!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More
{nat:>4d}   0  0  0.0 0.0 0.0   0   0   0   0   0        {atz}    0   0   0
!
{space_group}               <--Space group symbol
!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
""".format(number=number,
           name=name,
           nat=len(scatt),
           atz=0.0,
           space_group=symm.space_group_info() )
  for atom in scatt:
    tmp += """\
{lbl:s}     {stype:s}      {x:f}  {y:f}  {z:f}  {biso:7.5f}  {occ:7.5f}   0   0   0    0
             ##_ap*_## ##_ap*_## ##_ap*_##     0.00      0.00\n""".format(
                                      lbl=atom.label,
                                      stype=atom.element_symbol(),
                                      x=atom.site[0],
                                      y=atom.site[1],
                                      z=atom.site[2],
                                      biso=atom.b_iso(),
                                      occ=atom.occupancy)
  a, b, c, alpha, beta, gamma = symm.unit_cell().parameters()
  tmp += """\
!-------> Profile Parameters for Pattern #  1
!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model
  {scale:f}       0.10033   0.03837   0.00000   0.00000   0.00000      0
##_scf_##    ##_shp1.{n}_## ##_shp2.{n}_##   0.000    0.000     0.000
!       U         V          W           X          Y        GauSiz   LorSiz Size-Model
   1.706310  -1.179299   0.405540   0.000000   0.000000   0.000000   0.000000    0
 ##_prfu.{n}_## ##_prfv.{n}_##  ##_prfw.{n}_##      0.000      0.000      0.000      0.000
!     a          b         c        alpha      beta       gamma      #Cell Info
   {a:f}  {b:f}  {c:f}  {alpha:f}  {beta:f}   {gamma:f}\n""".format(
                                    scale=0.01/scale_down,
                                    a=a, b=b, c=c,
                                    alpha=alpha, beta=beta, gamma=gamma, n=number)
  # make contraints from symmetry
  def __constraints_from_symm(xtal_system):
    xtal_system = xtal_system.lower()
    if xtal_system == "triclinic":
      return "  ##_lp1.{0}_## ##_lp2.{0}_## ##_lp3.{0}_## ##_lp4.{0}_## ##_lp5.{0}_## ##_lp6.{0}_##".format(number)
    elif xtal_system == "monoclinic":
      return "  ##_lp1.{0}_## ##_lp2.{0}_## ##_lp3.{0}_## ##_lp4.{0}_##     0.00 ##_lp5.{0}_##".format(number)
    elif xtal_system == "orthorhombic":
      return "  ##_lp1.{0}_## ##_lp2.{0}_## ##_lp3.{0}_##     0.00     0.00     0.00".format(number)
    elif xtal_system in ["trigonal", "hexagonal", "tetragonal"]:
      return "  ##_lp1.{0}_## ##_lp1.{0}_## ##_lp2.{0}_##     0.00     0.00     0.00".format(number)
    elif xtal_system == "cubic":
      return "  ##_lp1.{0}_## ##_lp1.{0}_## ##_lp1.{0}_##     0.00     0.00     0.00".format(number)
    else:
      return "     0.00     0.00     0.00     0.00     0.00     0.00"
  tmp += __constraints_from_symm(symm.space_group().crystal_system()) + "\n"
  tmp += """\
!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4
  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
     0.00     0.00     0.00     0.00     0.00     0.00
"""
  return tmp


def _pcr_skelleton(phases,
                  title="unnamed",
                  jobtype=0,
                  nprof=0,
                  nbckgd=0,
                  wavelength=wavelengths.characteristic("CU").as_angstrom(),
                  scale_down=1.0):
  """Create a pcr skelleton with placeholder strings for different refinement
  options.

  (for internal use)

  :param phases: a list/tuple/set of structures
  :type phases: list(cctbx.xray_structure)
  :param title: a title to be used in the pcr file
  :type title: string
  :param jobtype: the jobtype to be used
  :type jobtype: integer
  :param nprof: the default peak profile to be used
  :type nprof: integer
  :param nbckgd: the type of background to be used
  :type nbckgd: integer
  :param wavelength: the wavelength to be used in Angstroms
  :type wavelength: float
  :param scale_down: factor to divide intensities by (to avoid overflows)
  :type scale_down: float

  :returns: the pcr skelleton as a string
  :rtype: string
  """
  nphase = len(phases)

  ret = "COMM " + title + "\n"
  ret += """\
!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut
  {0}   {1}   {2}   {3}   0   0   0   {dum}   0   0   0   0   0   0   0   0   0   0   1
""".format(jobtype,nprof,nphase,nbckgd,dum=1 if jobtype == 2 else 0)
  ret += """\
!
!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 Syo Prf Ins Rpa Sym Hkl Fou Sho Ana
  2   0   1   1   1   0   0   0   1   1  {filetype}   0   1   4   2   0   0
!
! lambda1 Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz ->Patt# 1
 {xlambda} {xlambda}  1.0000   60.000  5.0000  0.0000  0.0000   50.00    0.0000
!
!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0
 20  0.01  1.00  1.00  1.00  1.00     10.0000   0.100000   100.0000   0.000   0.000
!
!
  #__npar__#    !Number of refined parameters
""".format(filetype="0", xlambda=wavelength)
  ret += """\
!
!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1
  0.00000    0.0  0.00000    0.0  0.00000    0.0 0.000000    0.00   0
!   Background coefficients/codes  for Pattern#  1  (Polynomial of 6th degree)
       0.000       0.000       0.000       0.000       0.000       0.000
        0.00        0.00        0.00        0.00        0.00        0.00
"""
  for i in range(nphase):
    ret += _make_phase_block(phases[i], number=i+1, name="", scale_down=scale_down)
  return ret



def _set_ref_flags(inputstring, freeparams=[]):
  """Parse a pcr skelleton to enable or disable the different refinement
  parameters. The placeholder string for a parameter starts with '##_' and ends
  with '_##'.

  (for internal use)

  Allowed values in the freeparams list are:
    * 'scale'   --> free all scale factors
    * 'lattice' --> free all lattice parameters
    * 'profile' --> free all profile parameters
  The parameters are freed according to their order inside the list.
  If an empty list is passed all parameters will be fixed.

  :param inputstring: a string containing placeholders for different refinement \
  variables
  :type inputstring: string
  :param freeparams: a list of parameter sets to free
  :type freeparams: list

  :returns: the parsed string
  :rtype: string
  """
  import re
  varcount = 1
  ret = inputstring
  param_to_var = dict()
  def __replace_match(string, match, var):
    mlen = len(match.group(0))
    replace = "{0:>d}1.00".format(var).rjust(mlen)
    return string[:match.start()] + replace + string[match.end():]
  if freeparams is not []:
    for param in freeparams:
      if param == "scale":
        # free all scale factors
        for m in re.finditer('##_scf_##', ret):
          ret = __replace_match(ret, m , varcount)
          varcount += 1
      elif param == "lattice":
        # free all lattice parameters
        for m in re.finditer('##_lp.*?_##', ret):
          if m.group(0) not in param_to_var:
            param_to_var[m.group(0)] = varcount
            varcount += 1
          ret = __replace_match(ret, m , param_to_var[m.group(0)])
      elif param == "profile":
        # free all profile parameters
        for m in re.finditer('##_prf*?_##', ret):
          if m.group(0) not in param_to_var:
            param_to_var[m.group(0)] = varcount
            varcount += 1
          ret = __replace_match(ret, m , param_to_var[m.group(0)])
      else:
        raise ValueError("unknown parameter type: '{0}'".format(param))
  # fix all still unhandled flags
  for m in re.finditer('##_.*?_##', ret):
    ret = ret[:m.start()] + "0.00".rjust(len(m.group(0))) + ret[m.end():]
  ret = ret.replace('#__npar__#', str(varcount-1))
  return ret

def write_pcr(s,
              phases,
              title="unnamed",
              jobtype=0,
              nprof=0,
              nbckgd=0,
              wavelength=wavelengths.characteristic("CU").as_angstrom(),
              I_obs=None,
              scale_down=1.0):
  """Write a pcr file to a file or IO buffer.

  :param s: a buffer/file to write to
  :type s: StringIO or filestream
  :param phases: a list/tuple/set of structures or a single crystal structucture
  :type phases: list(cctbx.xray_structure) or cctbx.xray_structure
  :param title: a title to be used in the pcr file
  :type title: string
  :param jobtype: the jobtype to be used (0=xray refine, 2=xray powdersim)
  :type jobtype: integer
  :param nprof: the default peak profile to be used
  :type nprof: integer
  :param nbckgd: the type of background to be used
  :type nbckgd: integer
  :param wavelength: the wavelength to be used in Angstroms
  :type wavelength: float
  :param I_obs: observed intensities
  :type I_obs: cctbx.miller
  :param scale_down: factor to divide intensities by (to avoid overflows)
  :type scale_down: float
  """
  # handle case of being called with only one phase
  if not isinstance(phases, (tuple, list, set)):
    phases = (phases, )

  # XXX Todo: handle case of given I_obs (= refine without a full profile)
  pcrfile = _pcr_skelleton(phases=phases,
                  title=title,
                  jobtype=jobtype,
                  nprof=nprof,
                  nbckgd=nbckgd,
                  scale_down=scale_down)
  if jobtype == 2:
    pcrfile = _set_ref_flags(pcrfile)
  else:
    pcrfile = _set_ref_flags(pcrfile, freeparams=["scale","lattice","profile"])
  print(pcrfile, file=s)


if __name__ == '__main__':
  # just a little test for debugging
  from StringIO import StringIO
  from cctbx import sgtbx
  from cctbx.development import random_structure
  xrs1 = random_structure.xray_structure(
        space_group_info=sgtbx.space_group_info(number=1),
        elements=["C"]*10,u_iso=0.005)
  xrs2 = random_structure.xray_structure(
        space_group_info=sgtbx.space_group_info(number=123),
        elements=["C"]*4,u_iso=0.005)
  #cif = xrs1.as_cif_simple()
  #print(cif)
  buffer = StringIO()
  write_pcr(buffer, (xrs1, xrs2), jobtype=2)
  print(buffer.getvalue())
  buffer.close
