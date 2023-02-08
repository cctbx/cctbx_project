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
from __future__ import absolute_import, division, print_function
import os
import re
import subprocess
from six.moves import cStringIO as StringIO
import sys
import time
import stat

from libtbx import easy_run
from libtbx.utils import Sorry
from six.moves import range

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
  "python_script"    : str,
  "parallel_nodes"   : int,
  'fake_queue'       : bool,
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

  #assert chunk_size==len(cmds) or chunk_n==1
  assert chunk_n>0
  assert chunk_size>0
  #if chunk_n!=1:
  #  chunk_size = (len(cmds)-1)//chunk_n+1
  #elif chunk_size!=1:
  #  chunk_n = len(cmds)%%chunk_size+1

  for i, cmd in enumerate(cmds):
    if only_i is None:
      print(i, cmd)
      continue
    else:
      if chunk_n!=1 or chunk_size!=len(cmds):
        if only_i!=i//chunk_size+1: continue
      else:
        if only_i!=i+1: continue
    print('Running', i+1)
    print('Command', cmd)
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

libtbx.python %s $SGE_TASK_ID %s %s >& %s.$SGE_TASK_ID.out

%s

exit

"""

test_run_script = """
import os, sys

def run(only_i=None):
  print("%s " % only_i)
  f=open("test_%s.output" % only_i, "w")
  f.write("Testing %s\\n" % only_i)
  f.close()

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
"""

qblock = """#! /bin/csh -q
#$ -j y
while (1)
  sleep 3600
end
"""

def run_sh(cmd):
  t0=time.time()
  rc = easy_run.call(cmd)
  print('done "%s" in %0.1f' % (cmd.strip(),
                                time.time()-t0))

def test_easy_qsub():
  def _clean():
    for filename in os.listdir(os.getcwd()):
      if filename.startswith("commands"):
        print('remove',filename)
        os.remove(filename)
  def _write_cmds(number_of_cmds):
    print("\n  test list of commands")
    outl = ""
    for i in range(number_of_cmds):
      outl += "echo %d\n" % (i+101)
    print(outl)
    f=open("commands.txt", "w")
    f.write(outl[:-1])
    f.close()
  def _wait():
    import time
    time.sleep(5)
    running=True
    while running:
      running=False
      for line in os.popen("qstat"):
        if line.find("commands")>-1:
          print('running',line)
          running=True
          break
  def _check(number_of_output_files,
             number_of_runs,
             ):
    runs=[]
    for filename in os.listdir(os.getcwd()):
      if filename.startswith("commands") and filename.endswith(".out"):
        number_of_output_files-=1
        f=open(filename, "r")
        lines=f.readlines()
        f.close()
        for line in lines:
          if line.find("Running")>-1:
            number_of_runs-=1
            tmp = line.split()
            assert tmp[1] not in runs
            runs.append(tmp[1])

    assert not number_of_output_files
    assert not number_of_runs

  print('#'*80)
  print('# testing easy_qsub')
  print('#'*80)
  cmd = 'libtbx.easy_qsub phenix_source="/net/cci/xp/phenix/phenix_env" commands=commands.txt'
  old_cmd = cmd
  for cmds, new_cmd, number_of_output_files, number_of_runs in [
      [10,"",10,10],
      [10," number_of_chunks=2", 2, 10],
      [11," number_of_chunks=2", 2, 11],
      [10," size_of_chunks=2",   5, 10],
      [11," size_of_chunks=2",   6, 11],
      [10," start=6",            5, 5],
      [10," end=5",              5, 5],
      [20," start=5 size_of_chunks=2", 6, 12],
      [20," end=5 size_of_chunks=2", 5, 10],
      ]:
    cmd = old_cmd + new_cmd
    _clean()
    _write_cmds(cmds)
    os.system(cmd)
    _wait()
    _check(number_of_output_files=number_of_output_files,
           number_of_runs=number_of_runs,
      )
  _clean()
  #
  sys.exit()

def process_args(args):
  kwds = {}
  for t in args:
    if t in ["--help", "-h"]:
      print("""
  Program to sumbit jobs easily to a SGE queue
    e.g.

    libtbx.easy_qsub phenix_source="/net/cci/xp/phenix/phenix_env" commands=commands.txt

  """)
      sys.exit()
    elif t=="--test":
      test_easy_qsub()
      return {}
    elif t=="--dry":
      kwds["dry_run"]=True
    for key in key_words:
      if t.find("%s=" % key)==0:
        kwds[key]=t.split("=")[1]
        kwds[key] = key_words[key](kwds[key])
        break
    else:
      if t not in ["--dry"]:
        print('\n  failed to process "%s"\n' % t)
        assert 0
  return kwds

def get_queue_machine_details():
  cmd = "qstat -f"
  ero = easy_run.fully_buffered(command = cmd)
  out = StringIO()
  ero.show_stdout(out = out)
  rc = {}
  for line in out.getvalue().split("\n"):
    if line.find("all.q@")>-1:
      tmp = line.split()
      rc.setdefault(tmp[0], [])
      rc[tmp[0]].append(int(tmp[2].split("/")[0]))
      rc[tmp[0]].append(int(tmp[2].split("/")[1]))
  return rc

def get_python_bin(source):
  import libtbx.load_env
  return libtbx.env.python_exe.sh_value()

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
        python_script=None,
        parallel_nodes=1,
        fake_queue=False,
        dry_run=False,
        ):
  help = """
  Submits a job to run commands on the queueing system using qsub. One master
  job is sent, which then deploys multiple sub jobs (Of the same job ID) to run
  each separate command.

  phenix_source: A path to the script used to initial the phenix environment
  where: A path to the directory in which to run the commands
  commands: A list of commands to run on the queue
  size_of_chunks: Number of commands to run within each queue job
  code: The task name used during submissions
  js: queue priority (default 0)
  qsub_cmd: The command to submit a job from a shell
  start: The start number for the task range
  end: The end number for the task range
  host_scratch_dir: If not None, will run the job within a local scratch
                    directory on each queue machine, then copy back all output
                    files to the directory set by where

  Returns the job id of the submission.

  Example of python script:
  run(
    where = path_to_where_you_want_to_run_your_jobs,
    phenix_source = "/net/chevy/raid1/nigel/build/setpaths.csh",
    commands = ["phenix.fetch_pdb 101m", "phenix.fetch_pdb 1s72"])

  For more help type:
    libtbx.easy_qsub --help
  """
  if not phenix_source:
    print('-'*80)
    if not os.environ.get("PHENIX", ""):
      print(help)
      return
    print("\n  Automatically setting phenix_source to current $PHENIX/phenix_env\n")
    phenix_source = "%s/phenix_env" % os.environ.get("PHENIX", "")

  if not (commands or python_script):
    print('-'*80)
    print("\n  Generating a test run script and queuing 10 times\n")
    f=open("easy_qsub_test_script.py", "w")
    f.write(test_run_script)
    f.close()
    commands = []
    for i in range(10):
      commands.append("libtbx.python %s %s" % (os.path.join(
                                               os.getcwd(),
                                               "easy_qsub_test_script.py"),
                                               i+100,
                                              )
        )


  print('-'*80)
  print('  Inputs')
  print('    phenix_source',phenix_source)
  if phenix_source.find("phenix_env")==-1 and phenix_source.find("setpath")==-1:
    print('  Need to supply file to source. e.g. phenix_env')
    return False
  if not os.path.exists(phenix_source):
    raise Sorry('source file for PHENIX environment not found "%s"' % phenix_source)
  print('    where',where)
  if isinstance(commands, type([])):
    if code is None: code = "easy_qsub"
    print('    commands',len(commands), end=' ')
    if len(commands)>1:
      print('similar to\n\n> %s\n' % (commands[0]))

  elif commands is not None and os.path.exists(commands):
    if code is None: code = commands[:10]

  elif python_script is not None and os.path.exists(python_script):
    if code is None: code = python_script.replace("../","")[:10]

  code = code.replace("/", "")

  print('    size_of_chunks',size_of_chunks)
  print('    number_of_chunks',number_of_chunks)
  if number_of_chunks==1:
    print('\n  Need to choose number_of_chunks>1')
    return
  print('-'*80)
  old_where = os.getcwd()
  if where is None:
    where = old_where

  assert phenix_source

  if isinstance(commands, type([])):
    lines = commands
  elif commands is not None and os.path.exists(commands):
    f=open(commands, "r")
    lines = f.readlines()
    f.close()

  elif python_script is not None and os.path.exists(python_script):
    if not end: raise Sorry('Need to supply "end=n" for python_script')
    phenix_python_bin = get_python_bin(phenix_source)
    assert phenix_python_bin
    lines = []
    for i in range(start, end+1):
      lines.append("%s %s %d" % (phenix_python_bin, python_script, i))

  number_of_jobs = len(lines)
  print('\n  Number of lines in command file',number_of_jobs)
  if number_of_chunks is not None:
    number_of_chunks = min(number_of_chunks, number_of_jobs)
  if number_of_chunks is None:
    if size_of_chunks==1:
      number_of_chunks=len(lines)
    else:
      number_of_chunks=number_of_jobs//size_of_chunks
      if number_of_jobs%size_of_chunks: number_of_chunks+=1
  else:
    size_of_chunks = number_of_jobs//number_of_chunks
    size_of_chunks = max(1, size_of_chunks)
    if number_of_chunks%size_of_chunks: size_of_chunks+=1

  number_of_jobs = number_of_chunks

  print('\n  Number of queue jobs',number_of_jobs)
  print('  Number of command in each queue job',size_of_chunks)

  if fake_queue:
    print('\nCreating fake queue with %s parallel nodes' % parallel_nodes)
    print(commands)
    f=open(commands, 'r')
    lines = f.readlines()
    f.close()
    from multiprocessing import Pool
    pool = Pool(parallel_nodes)
    pool.map(run_sh, lines)
    return None

  elif parallel_nodes>1:
    assert 0, 'parallel_nodes removed'
    def _cmp_gap(k1, k2):
      if details[k1][1]-details[k1][0]>details[k2][1]-details[k2][0]:
        return -1
      return 1
    #
    assert number_of_jobs==1, "Only one job can be run in parallel"
    details = get_queue_machine_details()
    keys = sorted(details, cmp=_cmp_gap)
    for queue_name in keys:
      if parallel_nodes<=details[queue_name][1]-details[queue_name][0]:
        # submit
        f=open("qblock.csh", "w")
        f.write(qblock)
        f.close()
        os.chmod("qblock.csh", stat.S_IREAD|stat.S_IWRITE|stat.S_IXUSR)
        qsub_cmd += " -q %s" % queue_name
        cmd = "%s -t 1-%d qblock.csh" % (qsub_cmd, parallel_nodes-1)
        print("  Blocking %d slots on %s" % (parallel_nodes-1, queue_name))
        print("    Need to remove them manually")
        print("   ",cmd)
        easy_run.call(cmd)
        break
    else:
      raise Sorry("No queue machine found with %d available slots" % parallel_nodes)

  #print '\n  Changing to work directory : %s' % where
  os.chdir(where)

  if host_scratch_dir:
    print('\n  Setting up scratch directories on queue hosts')
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
  print("  Writing queue python script:\n    %s" % python_run_filename)
  f=open(python_run_filename, "w")
  f.write(script_file % outl)
  f.close()

  qsub_run_filename = os.path.join(os.getcwd(), qsub_run_filename)
  print("  Writing queue command script:\n    %s" % qsub_run_filename)
  f=open(qsub_run_filename, "w")
  f.write(run_file % (
    code,
    code,
    phenix_source,
    pre,
    python_run_filename,
    number_of_chunks,
    size_of_chunks,
    code,
    post,
    )
    )
  f.close()

  if end is not None:
    number_of_jobs = min(number_of_jobs, end)

  cmd = "%s -t %d-%d -js %d %s" % (
    qsub_cmd,
    start,
    number_of_jobs,
    js,
    qsub_run_filename,
    )
  print('\n  Queue command\n')
  print("> %s\n" % cmd)
  if dry_run:
    print('  Skipping run...')

  # Run the command, and then find the job id in the output
  ero = easy_run.fully_buffered(command = cmd)
  out = StringIO()
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

  This function returns when there are no more jobs matching any of the
  parameters. i.e. No jobs matching <user>, no jobs matching <job_id>, and no
  jobs matching <task_name>.

  Parameters
  ----------
  task_name : int or str or list of int or list of str or None, optional
  job_id : int or str or list of int or list of str or None, optional
  user : str, optional
  qstat_cmd : str, optional
  sleep_time : int, optional
  """

  # Turn job_id into an array of strings if it isn't already one
  if job_id:
    if hasattr(job_id, "__iter__"):
      job_id = [str(i).strip() for i in job_id]
    else:
      job_id = [str(job_id).strip()]

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
