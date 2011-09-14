import sys

def run(args):
  assert len(args) == 1
  s = open(args[0], "rb").read().replace("__SOURCEFORGE_LOGO_HERE__", """\
<a href="http://sourceforge.net/projects/cctbx"
><img src="http://sflogo.sourceforge.net/sflogo.php?group_id=24107&type=14"
width="150" height="40" border="0"
alt="Get Computational Crystallography Toolbox at SourceForge.net.
Fast, secure and Free Open Source software downloads" /></a>""")
  s = open(args[0], "wb").write(s)

if (__name__ == "__main__"):
  run(sys.argv[1:])
