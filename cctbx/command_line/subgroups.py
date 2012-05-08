def run(args):
  if (len(args) == 0):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage(
      "%s all|list-of-space-group-symbols-or-numbers"
        % libtbx.env.dispatcher_name)
  from cctbx.sgtbx import space_group_info
  from cctbx.sgtbx.subgroups import show
  if (args == ["all"]):
    for space_group_number in xrange(1,231):
      show(parent_group_info=space_group_info(space_group_number))
  else:
    for arg in args:
      show(parent_group_info=space_group_info(symbol=arg))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
