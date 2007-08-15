
#         1         2         3         4         5         6         7
#1234567890123456789012345678901234567890123456789012345678901234567890


def write_remark_3(info, out):
  pr = "REMARK   3 "
  print >> out,pr
  print >> out,pr+"REFINEMENT.                                                "
  print >> out,pr+" PROGRAM     : PHENIX (phenix.refine)                      "
  print >> out,pr+" AUTHORS     : Paul Adams, Pavel Afonine, Vicent Chen, Ian "
  print >> out,pr+"             : Davis, Kreshna Gopal, Ralf Grosse-Kunstleve,"
  print >> out,pr+"             : Li-Wei Hung, Robert Immormino, Tom Ioerger, "
  print >> out,pr+"             : Airlie McCoy, Erik McKee, Nigel Moriarty,   "
  print >> out,pr+"             : Reetal Pai, Randy Read, Jane Richardson,    "
  print >> out,pr+"             : David Richardson, Tod Romo, Jim             "
  print >> out,pr+"             : Sacchettini, Nicholas Sauter, Jacob Smith,  "
  print >> out,pr+"             : Laurent Storoni, Tom Terwilliger, Peter     "
  print >> out,pr+"             : Zwart                                       "
  print >> out,pr
  info.show_remark_3(out = out)
  print >> out,pr
