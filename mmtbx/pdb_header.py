import os, sys, libtbx
#         1         2         3         4         5         6         7
#1234567890123456789012345678901234567890123456789012345678901234567890

def version(f=None):
  if (f is None): f = sys.stdout
  version = os.environ.get("PHENIX_VERSION", None)
  if (version is None):
    tag_file = libtbx.env.under_dist("libtbx", "../TAG")
    if (os.path.isfile(tag_file)):
      try: version = open(tag_file).read().strip()
      except KeyboardInterrupt: raise
      except: pass
  release_tag = os.environ.get("PHENIX_RELEASE_TAG", None)
  cctbx_tag_file = libtbx.env.under_dist("libtbx", "../cctbx_bundle_TAG")
  cctbx_tag = None
  if (os.path.isfile(cctbx_tag_file)):
    try: cctbx_tag = open(cctbx_tag_file).read().strip()
    except KeyboardInterrupt: raise
    except: pass
  result = ""
  if(version is not None):
    result += version
  if (release_tag is not None):
    result = result + "-" +release_tag
  if(len(result.strip())==0): result = None
  return result


def write_remark_3(info, out):
  ver = version(f = out)
  if(ver is None):
    prog = " PROGRAM     : PHENIX (phenix.refine)                      "
  else:
    prog = " PROGRAM     : PHENIX (phenix.refine: %s)"%ver
  pr = "REMARK   3 "
  print >> out,pr
  print >> out,pr+"REFINEMENT.                                                "
  print >> out,pr+prog
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
