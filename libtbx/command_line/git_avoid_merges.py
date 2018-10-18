from __future__ import absolute_import, division, print_function
from libtbx.auto_build.bootstrap import Toolbox
import libtbx.easy_run
import libtbx.load_env
import os
import sys

what_is_this = '''

The command libtbx.git_avoid_merges identifies all module directories that are
git repositories (ie. it checks all paths listed by libtbx.show_dist_paths),
and changes the git settings of those repositories to 'rebase' by default
instead of creating merge commits. The global git settings are also changed
so any new git repositories are also set to 'rebase' by default.
Any existing settings are left as they are.

What does 'rebase' vs 'merge' mean?
A 'git pull' (without rebase setting or '--rebase' option) will silently
introduce a new commit, a so-called merge commit, if you and someone else have
committed anything since the last time you pulled. Invariably this will create
a lot of merge commits, going quite far back in history as well. Merge commits
will generally have a useless commit message, clutter the project history and
make it difficult to see what was going on. With our repositories the merge
commits will also cause pointless notification mail spam.

A 'rebase' will take your commits off the current tree and tack them back on to
the most recent commit. This results in a linear history.

Now there may be situations where you want to override this default policy.
You can do that by specifying --no-rebase on the command line, eg.:
   git pull --no-rebase

'''

def find_all_git_modules():
  for path in (abs(p) for p in libtbx.env.repository_paths):
    for entry in os.listdir(path):
      if not entry.startswith('.'):
        modulepath = os.path.join(path, entry)
        configfile = os.path.join(modulepath, '.git', 'config')
        if os.path.exists(configfile):
          yield (modulepath, configfile)

def mangle_git_repository(module, config):
  t = Toolbox()
  print("Git repository found:", module)
  t.set_git_repository_config_to_rebase(config)

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
    print("Updating global git settings:")
    for attribute, value in expected.items():
      print("  %s = %s" % (attribute, value))
      libtbx.easy_run.fully_buffered("git config --global \"%s\" \"%s\"" % (attribute, value))

def run():
  for module, config in find_all_git_modules():
    mangle_git_repository(module, config)
  set_git_defaults_to_rebase()
  print("\nAll done.")

if __name__ == "__main__":
  if len(sys.argv) != 1:
    print(what_is_this)
  else:
    run()
