import sys, os

def run(args):
  for path in args:
    problem = None
    if (not os.path.isfile(path) or os.path.islink(path)):
      problem = "not a regular file"
    else:
      try:
        file_content = open(path, "rb").read()
      except Exception:
        problem = "no read access"
      else:
        if (not os.access(path, os.W_OK)):
          problem = "no write access"
    if (problem is not None):
      print "%s: %s -> no action" % (path, problem)
    else:
      n_cr = file_content.count("\r")
      n_lf = file_content.count("\n")
      n_crlf = file_content.count("\r\n")
      action = "unknown -> no action"
      unix_content = None
      if (n_crlf > 0 and n_crlf == n_cr):
        action = "dos -> unix"
        unix_content = file_content.replace("\r\n", "\n")
        if (ord(unix_content[-1]) == 26):
          unix_content = unix_content[:-1]
      elif (n_cr > 0 and n_lf == 0):
        action = "mac -> unix"
        unix_content = file_content.replace("\r", "\n")
      elif (n_lf > 0 and n_cr == 0):
        action = "unix -> no action"
      print "%s: %s" % (path, action)
      if (unix_content is not None):
        if (unix_content[-1] != "\n"):
          unix_content += "\n"
        try:
          open(path, "wb").write(unix_content)
        except Exception:
          print >> sys.stdout, "FATAL ERROR: Cannot write file:", path
          path_copy = path + "_copy"
          print >> sys.stdout, "Saving copy of old content as file:", path_copy
          open(path_copy, "wb").write(file_content)
          sys.exit(1)

if (__name__ == "__main__"):
  run(sys.argv[1:])
