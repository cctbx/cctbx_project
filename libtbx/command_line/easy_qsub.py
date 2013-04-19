"""
libtbx.qsub(
 where=path_to_where_you_want_to_run_your_jobs,
 source="/net/chevy/raid1/afonine/build/setpaths.csh",
 commands=[list of commands])

where [list of commands]) is

["python run.py model1.pdb data1.mtz",
 "python run.py model2.pdb data2.mtz",
...,
"python run.py modelN.pdb dataN.mtz" ]
"""
from __future__ import division
import getpass
import os
import re
import subprocess
import StringIO
import sys
import time

from libtbx import easy_run
from libtbx.utils import Sorry

key_words = {
  "phenix_source"    : str,
  "where"            : str,
  "commands"         : str,
  "number_of_chunks" : int,
  "size_of_chunks"   : int,
  "code"             : str,
  "js"               : int,
  "qsub_cmd"         : str,
  "start"            : int,
  "end"              : int,
  "host_scratch_dir" : str,
  }

script_file = """
import os, sys
from libtbx import easy_run

cmds = [
%s
]

def run(only_i=None,
        chunk_n=1,
        chunk_size=len(cmds),
        ):
  try: only_i = int(only_i)
  except ValueError: only_i=None

  try: chunk_n = int(chunk_n)
  except ValueError: chunk_n=1
  try: chunk_size = int(chunk_size)
  except ValueError: chunk_size=len(cmds)

  assert chunk_size==len(cmds) or chunk_n==1
  assert chunk_n>0
  assert chunk_size>0
  if chunk_n!=1:
    chunk_size = (len(cmds)-1)//chunk_n+1
  elif chunk_size!=1:
    chunk_n = len(cmds)%%chunk_size+1

  for i, cmd in enumerate(cmds):
    if only_i is None:
      print i, cmd
      continue
    else:
      if chunk_n!=1 or chunk_size!=len(cmds):
        if only_i!=i//chunk_size+1: continue
      else:
        if only_i!=i+1: continue
    print 'Running', i+1
    print 'Command', cmd
    easy_run.call(cmd)

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
"""

host_scratch_dir_pre = """
set work_dir=$PWD

set host_scratch_dir=/net/$host/scratch1/$user
if (! -d $host_scratch_dir) then
  echo "================================"
  echo " Creating host scratch directory $host_scratch_dir"
  echo "================================"
  mkdir $host_scratch_dir
endif

set host_scratch_dir=/net/$host/scratch1/$user/%s
if (! -d $host_scratch_dir) then
  echo "================================"
  echo " Creating host scratch directory $host_scratch_dir"
  echo "================================"
  mkdir $host_scratch_dir
endif

set host_scratch_dir=/net/$host/scratch1/$user/%s/$SGE_TASK_ID
if (! -d $host_scratch_dir) then
  echo "================================"
  echo " Creating host scratch directory $host_scratch_dir"
  echo "================================"
  mkdir $host_scratch_dir
endif

echo "=============="
echo " Changing into $host_scratch_dir"
echo "=============="
cd $host_scratch_dir
"""
host_scratch_dir_post = """
echo "==================="
echo " Move files back to $work_dir"
echo "==================="
mv * $work_dir
"""

run_file = """#! /bin/csh -q
#$ -cwd
#$ -o %s_queue.output -j y -N %s
limit datasize 2000000

source %s

%s

phenix.python %s $SGE_TASK_ID $SGE_TASK_LAST >& %s.$SGE_TASK_ID.out

%s

exit

"""

test_run_script = """
import os, sys

def run(only_i=None):
  print "%s " % only_i
  f=file("test_%s.output" % only_i, "wb")
  f.write("Testing %s\\n" % only_i)
  f.close()

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
"""

def process_args(args):
  kwds = {}
  for t in args:
    for key in key_words:
      if t.find("%s=" % key)==0:
        kwds[key]=t.split("=")[1]
        kwds[key] = key_words[key](kwds[key])
        break
    else:
      print '\n  failed to process "%s"\n' % t
      assert 0
  return kwds

def run(phenix_source=None,
        where=None,
        commands=None,
        size_of_chunks=1,
        number_of_chunks=None,
        code=None,
        js=0,
        qsub_cmd="qsub",
        start=1,
        end=None,
        host_scratch_dir=None,
        ):
  """
  Submits a job to run commands on the queueing system using qsub. One master
  job is sent, which then deploys multiple sub jobs (Of the same job ID) to run
  each separate command.

  phenix_source: A path to the script used to initial the phenix environment
  where: A path to the directory in which to run the commands
  commands: A list of commands to run on the queue
  size_of_chunks: Number of commands to run within each queue job
  code: The task name used during submissions
  js: ???
  qsub_cmd: The command to submit a job from a shell
  start: The start number for the task range
  end: The end number for the task range
  host_scratch_dir: If not None, will run the job within a local scratch
                    directory on each queue machine, then copy back all output
                    files to the directory set by where

  Returns the job id of the submission.

  Ex:
  run(
    where = path_to_where_you_want_to_run_your_jobs,
    source = "/net/chevy/raid1/afonine/build/setpaths.csh",
    commands = ["phenix.fetch_pdb abcd", "phenix.fetch_pdb efgh"])
  """
  if not phenix_source:
    print '-'*80
    print "\n  Automatically setting phenix_source to current $PHENIX/phenix_env\n"
    phenix_source = "%s/phenix_env" % os.environ.get("PHENIX", "")

  if not commands:
    print '-'*80
    print "\n  Generating a test run script and queuing 10 times\n"
    f=file("easy_qsub_test_script.py", "wb")
    f.write(test_run_script)
    f.close()
    commands = []
    for i in range(10):
      commands.append("phenix.python %s %s" % (os.path.join(
                                               os.getcwd(),
                                               "easy_qsub_test_script.py"),
                                               i+100,
                                              )
        )


  print '-'*80
  print '  Inputs'
  print '    phenix_source',phenix_source
  if phenix_source.find("phenix_env")==-1 and phenix_source.find("setpath")==-1:
    print '  Need to supply file to source. e.g. phenix_env'
    return False
  if not os.path.exists(phenix_source):
    raise Sorry('source file for PHENIX environment not found "%s"' % phenix_source)
  print '    where',where
  if type(commands)==type([]):
#  if commands is None:
#    print '  Need to supply a list of commands, either by file or python list'
#    return False
#  elif type(commands)==type([]):
    if code is None: code = "easy_qsub"
    print '    commands',len(commands),
    if len(commands)>1:
      print 'similar to\n\n> %s\n' % (commands[0])
  else:
    print '    commands',commands
    assert os.path.exists(commands)
    if code is None: code = commands[:8]
  print '    size_of_chunks',size_of_chunks
  print '    number_of_chunks',number_of_chunks
  if number_of_chunks==1:
    print '\n  Need to choose number_of_chunks>1'
    return
  print '-'*80
  old_where = os.getcwd()
  if where is None:
    where = old_where

  assert phenix_source

  if type(commands)==type([]):
    lines = commands
  else:
    f=file(commands, "rb")
    lines = f.readlines()
    f.close()
  number_of_jobs = len(lines)
  print '\n  Number of lines in command file',number_of_jobs
  if number_of_chunks is None:
    if size_of_chunks==1:
      number_of_chunks=len(lines)
    else:
      number_of_chunks=number_of_jobs//size_of_chunks
      if number_of_jobs%size_of_chunks: number_of_chunks+=1
  else:
    size_of_chunks = number_of_jobs//number_of_chunks
    if number_of_jobs%size_of_chunks: size_of_chunks+=1
  number_of_jobs = number_of_chunks

  print '\n  Number of queue jobs',number_of_jobs
  print '  Number of command in each queue job',size_of_chunks

  #print '\n  Changing to work directory : %s' % where
  os.chdir(where)

  if host_scratch_dir:
    print '\n  Setting up scratch directories on queue hosts'
    pre = host_scratch_dir_pre % (host_scratch_dir, host_scratch_dir)
    post = host_scratch_dir_post
  else:
    pre = ""
    post = ""

  python_run_filename = "easy_qsub_python_script.py"
  qsub_run_filename = "easy_qsub_qsub_script.sh"

  outl = ""
  for line in lines:
    if line[-1]=="\n": line = line[:-1]
    outl += "  '''%s ''',\n" % line
  python_run_filename = os.path.join(os.getcwd(), python_run_filename)
  print "  Writing queue python script:\n    %s" % python_run_filename
  f=file(python_run_filename, "wb")
  f.write(script_file % outl)
  f.close()

  qsub_run_filename = os.path.join(os.getcwd(), qsub_run_filename)
  print "  Writing queue command script:\n    %s" % qsub_run_filename
  f=file(qsub_run_filename, "wb")
  f.write(run_file % (
    code,
    code,
    phenix_source,
    pre,
    python_run_filename,
    code,
    post,
    )
    )
  f.close()

  if end is None:
    end = number_of_jobs

  cmd = "%s -t %d-%d -js %d %s" % (
    qsub_cmd,
    start,
    end,
    js,
    qsub_run_filename,
    )
  print '\n  Queue command\n'
  print "> %s\n" % cmd

  # Run the command, and then find the job id in the output
  ero = easy_run.fully_buffered(command = cmd)
  out = StringIO.StringIO()
  ero.show_stdout(out = out)
  # Line looks like:
  # "Your job-array 5436256.1-122:1 ("easy_qsub") has been submitted"
  pattern = re.compile(
      r"Your job-array (\d+)\..* \(\".*\"\) has been submitted")

  for line in out.getvalue().split("\n"):
    if not line: continue

    match = pattern.search(line)
    if match:
      job_id = int(match.group(1))
      break
  else:
    print("Unable to determine job ID")
    job_id = None

  os.chdir(old_where)

  return job_id

def wait(task_name = "easy_qsub",
         job_id = None,
         user = None,
         qstat_cmd = "qstat",
         sleep_time = 10):
  """
  Repeatedly checks the queue and waits until all jobs for a user have
  completed. The caller may specify further information, such as the job name,
  id, and user if jobs of another user are to be monitored.

  job_id and task_name may be ints, strings, or lists of ints or strings.
  """

  if user is None:
    user = getpass.getuser()

  # Turn job_id into an array of strings if it isn't already one
  if job_id:
    if hasattr(job_id, "__iter__"):
      job_id = [str(i).strip() for i in job_id]
    else:
      job_id = [str(task_name).strip()]

  # Do the same for task_name
  if task_name:
    if hasattr(task_name, "__iter__"):
      task_name = [str(i).strip() for i in task_name]
    else:
      task_name = [str(task_name).strip()]

  # Monitor the queue until we don't see any more tasks for the user
  # or with a given job id / task name
  while True:
    p = subprocess.Popen([qstat_cmd], stdout = subprocess.PIPE)
    stdout = p.communicate()[0]
    for line in stdout.split("\n")[2:]:
      if not line: continue

      # Check for a job with the same ID
      if job_id and line[:7].strip() in job_id:
        # Break the for loop
        break

      if line[27:39].strip() == user:
        if task_name is not None:
          if line[16:26].strip() in task_name:
            # Break the for loop
            break
        else:
          # Break the for loop
          break
    else:
      # No lines matched, so break the while loop
      break

    time.sleep(sleep_time)

if __name__=="__main__":
  args=sys.argv[1:]
  kwds = process_args(args)
  run(**kwds)
