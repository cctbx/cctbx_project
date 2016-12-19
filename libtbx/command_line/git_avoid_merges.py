from __future__ import division
import libtbx.easy_run
import libtbx.load_env
import os
import re

def find_all_git_modules():
  for path in (abs(p) for p in libtbx.env.repository_paths):
    for entry in os.listdir(path):
      if not entry.startswith('.'):
        modulepath = os.path.join(path, entry)
        configfile = os.path.join(modulepath, '.git', 'config')
        if os.path.exists(configfile):
          yield (modulepath, configfile)

def set_all_git_module_branches_to_rebase():
  for module, config in find_all_git_modules():
    print "Git repository found:", module
    with open(config, 'r') as fh:
      cfg = fh.readlines()

    branch, remote, rebase = False, False, False
    insertions = []
    for n, line in enumerate(cfg):
      if line.startswith('['):
        if branch and remote and not rebase:
          insertions.insert(0, (n, branch))
        if line.startswith('[branch'):
          branch = line.split('"')[1]
        else:
          branch = False
        remote, rebase = False, False
      if re.match('remote\s*=\s*', line.strip()):
        remote = True
      if re.match('rebase\s*=\s*', line.strip()):
        rebase = True
    if branch and remote and not rebase:
      insertions.insert(0, (n + 1, branch))
    for n, branch in insertions:
      print "  setting branch %s to rebase" % branch
      cfg.insert(n, '\trebase = true\n')
    with open(config, 'w') as fh:
      fh.write("".join(cfg))

def set_git_defaults_to_rebase():
  result = libtbx.easy_run.fully_buffered("git config --global --list")
  if result.return_code:
    '''Could not run git executable. Bail out.'''
    return

  expected = {
    'branch.autosetuprebase': 'always',
    'pull.rebase': 'true',
  }
  for line in result.stdout_lines:
    if '=' in line:
      name, value = line.split('=', 1)
      if name in expected:
        del(expected[name])
  if expected:
    print "Updating global git settings:"
    for attribute, value in expected.iteritems():
      print "  %s = %s" % (attribute, value)
      libtbx.easy_run.fully_buffered("git config --global \"%s\" \"%s\"" % (attribute, value))

def run():
  set_all_git_module_branches_to_rebase()
  set_git_defaults_to_rebase()
  print "\nAll done."

if __name__ == "__main__":
  run()
