if (__name__ == "__main__"):
  from cctbx.omz.cod_select_and_pickle import run
  import libtbx.load_env
  import sys
  run(args=sys.argv[1:], command_name=libtbx.env.dispatcher_name)
