
#         1         2         3         4         5         6         7
#1234567890123456789012345678901234567890123456789012345678901234567890


def write_remark_3(info, out):
  pr = "REMARK   3 "
  print >> out,pr
  print >> out,pr+"REFINEMENT.                                                "
  print >> out,pr+" PROGRAM     : PHENIX (phenix.refine)                      "
  print >> out,pr+" AUTHORS     : Paul Adams, Pavel Afonine, Vincent Chen, Ian"
  print >> out,pr+"             : Davis, Kreshna Gopal, Ralf Grosse-Kunstleve,"
  print >> out,pr+"             : Jeffrey Headd, Li-Wei Hung, Robert          "
  print >> out,pr+"             : Immormino, Tom Ioerger, Airlie McCoy, Erik  "
  print >> out,pr+"             : McKee, Nigel Moriarty, Reetal Pai,  Randy   "
  print >> out,pr+"             : Read, Jane Richardson, David Richardson, Tod"
  print >> out,pr+"             : Romo, Jim Sacchettini, Nicholas Sauter,     "
  print >> out,pr+"             : Jacob Smith, Laurent Storoni, Tom           "
  print >> out,pr+"             : Terwilliger, Peter Zwart                    "
  print >> out,pr
  info.show_remark_3(out = out)
  print >> out,pr
