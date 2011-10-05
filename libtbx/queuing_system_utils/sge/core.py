#
# found here (09-29-10):
#   http://pastebin.com/4x9YWB3r
#
# Original copyright notice:
#
# sge.py
#
# Creative Commons Attribution License
# http://creativecommons.org/licenses/by/2.5/
#
# Trevor Strohman
#    First release: December 2005
#    Second version:  11 January 2006
# Jacob Biesinger
#    This version:  September 2010
#
# Bug finders: Fernando Diaz
#

"""Submits jobs to Grid Engine.
Handles job dependencies, output redirection, and job script creation.
Also mimics multiprocessing.Pool mapping functions and a reduce function.
If SGE is not available, uses multiprocessing module.

USAGE:
>>> # SGE Job and JobGroup, with dependencies
>>> firstJob = sge.Job('test_single', 'echo Hello World')
>>> secondJob = sge.JobGroup('test_array', 'echo Hello from $arg1 and $arg2', arguments=dict(arg1=map(str, range(5)), arg2=['foo', 'bar']))
>>> secondJob.addDependency(firstJob)
>>> sge.build_submission('.', [firstJob, secondJob])  # order in list doesn't matter
# output only if SGE not installed-- otherwise it's sent to stdout/$job_name_$ID
Hello World
Hello from 0 and foo
Hello from 1 and foo
Hello from 2 and foo
Hello from 4 and foo
Hello from 3 and foo
Hello from 0 and bar
Hello from 1 and bar
Hello from 2 and bar
Hello from 3 and bar
Hello from 4 and bar
>>> # SGE map and reduce
>>> sgepool = sge.SGEPool()
>>> from math import sqrt
>>> results = sgepool.map(sqrt, [i**2 for i in range(10)])  # runs locally or on SGE and retrieves result
>>> results
[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
>>> from operator import add
>>> sgepool.reduce_unordered(add, results)
45.0
>>> sum(results)
45.0
"""


import libtbx.load_env

setpaths = libtbx.env.under_build("setpaths.sh")
headerSGEOnly = "source %s\n" % setpaths
#headerAllJobs = None
headerAllJobs = """\
python_cmd=libtbx.python;
code_dir=.;
output_dir=.;
export PYTHONPATH=$PYTHONPATH:$code_dir;
"""

import subprocess
import os
import os.path
import time
from multiprocessing import Pool, cpu_count
import cPickle as pickle
import tempfile
import itertools

class Job:
  def __init__( self, name, command, queue=None ):
    self.name = name
    self.queue = queue
    self.command = command
    self.script = command
    self.dependencies = []
    self.submitted = 0

  def addDependency( self, job ):
    self.dependencies.append( job )

  def wait( self ):
    if not check_for_qsub():
      # dont know status of local jobs
      # TODO could have a reference to local pool and join() here or have
      # a job-specific async_result
      return
    finished = 0
    interval = 3
    while not finished:
      finished = self.isFinished()
      if not finished:
        time.sleep(interval)
        interval = min( interval + 5, 60 )

  def isFinished(self):
    if not check_for_qsub():
      return True
    return subprocess.call(args="qstat -j %s > /dev/null 2>&1" % (self.name),
      shell=True)

class JobGroup:
  def __init__( self, name, command, queue=None, arguments={} ):
    self.name = name
    self.queue = queue
    self.command = command
    self.dependencies = []
    self.submitted = 0
    self.arguments = arguments
    self.generateScript()

  def generateScript( self ):
    self.script = ""
    # total number of jobs in this group
    total = 1

    # for now, SGE_TASK_ID becomes TASK_ID, but we base it at zero
    self.script += """let "TASK_ID=$SGE_TASK_ID - 1"\n"""

    # build the array definitions
    for key in self.arguments.keys():
      values = self.arguments[key]
      line = ("%s_ARRAY=( " % (key))
      for value in values:
        line += "\'"
        line += value
        line += "\' "
      line += " )\n"
      self.script += line
      total *= len(values)
    self.script += "\n"

    # now, build the decoding logic in the script
    for key in self.arguments.keys():
      count = len(self.arguments[key])
      self.script += """let "%s_INDEX=$TASK_ID %% %d"\n""" % ( key, count )
      self.script += """%s=${%s_ARRAY[$%s_INDEX]}\n""" % ( key, key, key )
      self.script += """let "TASK_ID=$TASK_ID / %d"\n""" % ( count )

    # now, run the job
    self.script += "\n"
    self.script += self.command
    self.script += "\n"

    # set the number of tasks in this group
    self.tasks = total

  def addDependency( self, job ):
    self.dependencies.append( job )

  def wait( self ):
    if not check_for_qsub():
      # dont know status of local jobs
      # TODO could have a reference to local pool and join() here or have
      # a job-specific async_result
      return

    finished = 0
    interval = 3
    while not finished:
      finished = self.isFinished()
      if not finished:
        time.sleep(interval)
        interval = min( 2 * interval, 60 )

  def isFinished(self):
    if not check_for_qsub():
      return True
    return subprocess.call(args="qstat -j %s > /dev/null 2>&1" % (self.name),
      shell=True)


class SGEPool:
  '''Mimics multiprocessing.Pool in an SGE context, allowing python functions to
  be run on SGE without leaving python or having to create separate command-line calls.
  If SGE (qsub) is unavailable, then jobs are run locally.'''

  def __init__(self, initializer=None, initargs=None, use_grid_engine=True):
    self.initializer = initializer
    self.initargs = initargs
    self.use_grid_engine = use_grid_engine
    if not check_for_qsub():
      # grid engine not available
      self.use_grid_engine = False

  def map(self, func, iterable, chunksize=1):
    '''An SGE equivalent of the map() built-in function (only one iterable argument), blocking until the result is ready.
    iterable can be divided into chunksize-d pieces
    '''
    if not self.use_grid_engine:  #
      workerPool = Pool(initializer=self.initializer, initargs=self.initargs)
      return workerPool.map(func, iterable, chunksize)
    else:
      return list(self.imap(func, iterable, chunksize))

  def imap(self, func, iterable, chunksize=1):
    '''SGE equivalent of map() built-in function (only one iterable argument).
    Note: All jobs are submitted simultaneously, but the results are returned lazily.
    '''
    if not self.use_grid_engine:
      workerPool = Pool(initializer=self.initializer, initargs=self.initargs)
      for val in workerPool.imap(func, iterable, chunksize):
        yield val
    else:
      iterable = iter(iterable)
      allJobs = self._submit_jobs(func, iterable, 'map', chunksize)
      #print len(allJobs)
      for job in allJobs:
        job.wait()  # wait for the job to finish
        for data in self._getData(job) : #.outputFile):
          yield data
        os.remove(job.inputFile)    # BUG: these files aren't removed if there is an exception raised
        os.remove(job.outputFile)

  def imap_unordered(self, func, iterable, chunksize=1):
    '''Same as SGEPool.imap, except that the results are unordered.
    Rather than blocking to ensure the correct order, all jobs are polled and
    results are returned as soon as they are done.
    '''
    if not self.use_grid_engine:
      workerPool = Pool(initializer=self.initializer, initargs=self.initargs)
      for val in workerPool.imap_unordered(func, iterable, chunksize):
        yield val
    iterable = iter(iterable)
    allJobs = self._submit_jobs(func, iterable, 'map', chunksize)
    interval = 3
    while len(allJobs) > 0:
      doneJobs = []
      for job in allJobs:
        if job.isFinished():
          doneJobs.append(job)
          for data in self._getData(job) : #.outputFile):
            yield data
          os.remove(job.inputFile)    # BUG: these files aren't removed if there is an exception raised
          os.remove(job.outputFile)
      for job in doneJobs:
        allJobs.remove(job)
      if len(doneJobs) == 0:
        # no jobs are done yet-- wait for a while for them to finish
        time.sleep(interval)
        interval = min( 2 * interval, 60 )

  def reduce_unordered(self, func, iterable, initializer=None, chunksize=2):
    '''Use SGE to apply a function of *two arguments* to the items of iterable, without regard to iterable's order.
    The final result is a single value, produced by applying the function to all the elements in iterable.
    For example, reduce_unordered(operator.add, [1,2,3,4,5]) would calculate /one/ of:
      ((((1+2)+3)+4)+5)   (1+(((2+3)+4)+5))  (1+(2+3)+(4+5))  etc.
    The SGE Jobs can handle more elements of the iterable by specifying chunksize >2.

    If the optional initializer is present, it is placed as the first element in each of the initial jobs.
    The initializer also serves as a default when the iterable is empty.
    If initializer is not given and iterable contains only one item, the first item is returned.

    For a chunksize of 3 with 10 items in the iterable, the job structure will look like:
    1  2  3  4  5  6  7  8  9  10
      \  |  /    \  |  /    \  |  /     \-----\
       Job 1      Job 2       Job 3        \ (wraps around to Job 4)
     \  |   /-------/        |
      Job 4          /------/
          \----\   /----/
            Job 5
    '''
    if chunksize != 0 and chunksize < 2:
      raise ValueError("chunksize in reduce_unordered must be >= 2, or 0 !")
    if not self.use_grid_engine:
      if initializer:
        return reduce(func, iterable, initializer)
      else:
        return reduce(func, iterable)
    iterable = iter(iterable)
    newJobs = []
    moreData = True
    # submit jobs of size chunksize until the iterable is exhausted,
    while moreData:
      if chunksize == 0:
        curData = list(iterable)
      else:
        curData = list(itertools.islice(iterable, chunksize))
      if len(curData) == 0:
        # iterable is exhausted of values
        moreData = False
      elif len(curData) == 1:
        # one more dataset to be included in next round of jobs
        newJobs.extend(curData)
      else:
        newJobs.extend(self._submit_jobs(func, curData, 'reduce', chunksize))
    # then create a layer of *new jobs* that depend on the previous set of jobs
    while len(newJobs) > 0:
      #print len(newJobs)
      # newJobs is of mixed datatype-- either SGE Job or original data.
      # If there are jobs, then the created jobs must depend on the old jobs
      curData = list(itertools.islice(newJobs, chunksize))
      #print curData
      if len(curData) == 1:
        # only one piece of data left-- completely reduced now
        if curData[0].__class__ == Job:
          curData[0].wait()  # wait for the last job in the hierarchy to finish
          reduceResult = pickle.load(open(curData[0].outputFile))
          newJobs.remove(curData[0])
        else:
          # no SGE jobs were run-- only one data point in iterable
          reduceResult = curData[0]
          newJobs.remove(curData[0])
      else:
        dependJobs = filter(lambda job: job.__class__ == Job, curData)
        addedJobs = self._submit_jobs(func, curData, 'reduce', chunksize, dependencies=dependJobs)
        for job in curData:
          newJobs.remove(job)
        newJobs.extend(addedJobs)
    return reduceResult

  def _submit_jobs(self, func, iterable, mode, chunksize=1, dependencies=[]):
    if mode not in ['map', 'reduce']:
      raise ValueError('mode must be one of "map" or "reduce"')
    assert (chunksize is not None)
    moreData = True
    allJobs = []
    iterable = iter(iterable)
    # for each partition
    while moreData:
      curChunks = itertools.islice(iterable, chunksize)
      # save the function and data to a new pickle file
      workFile = tempfile.NamedTemporaryFile(prefix='SGE_%s_'%mode, dir=os.path.abspath('.'), delete=False)
      pickle.dump(os.environ, workFile)   # dump all data about environment variables
      pickle.dump(self.initializer, workFile)
      if mode == 'map':
        pickle.dump(self.initargs, workFile)
      try:
        pickle.dump(func, workFile)     # dump mapping function
      except TypeError, e:
        # Failure pickling function-- probably a lambda function
        os.remove(workFile.name)
        raise e  # Can't pickle lambda or shell-defined functions
      curDataCount = 0
      for data in curChunks:
        curDataCount += 1
        pickle.dump(data, workFile)
      if curDataCount == 0:
        # no data was sliced from original dataset-- we are done
        moreData = False
        os.remove(workFile.name)
      else:
        # submit an SGE job on this dataset slice
        outfile = tempfile.NamedTemporaryFile(prefix='SGE_%s_results_'%mode, dir=os.path.abspath('.'), delete=False)
        cmd = 'libtbx.python "%s" --mode=%s "%s" "%s"' % (os.path.abspath(__file__), mode, workFile.name, outfile.name)
        curJob = Job(os.path.split(workFile.name)[1], cmd)
        for job in dependencies:
          curJob.addDependency(job)
        curJob.inputFile = workFile.name
        curJob.outputFile = outfile.name
        build_submission(os.path.abspath('.'), [curJob])
        allJobs.append(curJob)
    return allJobs

  def _getData (self, job) : #pickleFileName):
    '''Return all the data within a pickled file, checking for and re-raising any Exceptions found.'''
    with open(job.outputFile) as resultsFile:
      moreData = True
      while moreData:
        try:
          result = pickle.load(resultsFile)
        except EOFError:
          moreData = False
        else:
          if issubclass(type(result), Exception):
            raise result    # There was an error in mapping, which has been passed on to us
          yield result

def build_directories( directory ):
  subdirectories = [ "output", "stderr", "stdout", "jobs" ];
  #subdirectories = [ "stderr", "stdout", "jobs" ];
  directories = [ os.path.join( directory, subdir ) for subdir in subdirectories ]
  needed = filter( lambda x: not os.path.exists( x ), directories )
  map( os.mkdir, needed )

def build_job_scripts( directory, jobs, use_grid_engine ):
  for job in jobs:
    scriptPath = os.path.join( directory, "jobs", job.name )
    scriptFile = file( scriptPath, "w" )
    scriptFile.write( "#!/bin/bash\n" )
    scriptFile.write( "#$ -S /bin/bash\n" )
    if headerSGEOnly and use_grid_engine:
      scriptFile.write( headerSGEOnly )
    if headerAllJobs:
      scriptFile.write( headerAllJobs )
    #scriptFile.write( "cd %s \n" % os.path.join(directory,"output"))
    scriptFile.write( job.script + "\n" )
    scriptFile.close()
    os.chmod( scriptPath, 0755 )
    job.scriptPath = scriptPath

def extract_submittable_jobs( waiting ):
  """Return all jobs that aren't yet submitted, but have no dependencies that have
  not already been submitted."""
  submittable = []

  for job in waiting:
    unsatisfied = sum([(subjob.submitted==0) for subjob in job.dependencies])
    if unsatisfied == 0:
      submittable.append( job )

  return submittable

def submit_safe_jobs( directory, jobs ):
  """Submits a list of jobs to Grid Engine, using the directory structure provided to store
  output from stderr and stdout."""
  for job in jobs:
    job.out = os.path.join( directory, "stdout" )
    job.err = os.path.join( directory, "stderr" )

    args = " -N %s " % (job.name)
    args += " -o %s -e %s " % (job.out, job.err)

    if job.queue != None:
      args += "-q %s " % job.queue

    if isinstance( job, JobGroup ):
      args += "-t 1:%d " % ( job.tasks )

    if len(job.dependencies) > 0:
      args += "-hold_jid "
      for dep in job.dependencies:
        args += dep.name + ","
      args = args[:-1]

    qsubcmd = ("qsub %s %s" % (args, job.scriptPath))
    #print qsubcmd
    subprocess.call(args=qsubcmd, shell=True)
    job.submitted = 1

def run_safe_jobs( directory, jobs ):
  """In the event that Grid Engine is not installed, this program will
  run jobs serially on the local host."""
  pool = Pool(processes=max(cpu_count(), 1))
  errorCodes = []
  for job in jobs:
    job.out = os.path.join( directory, "stdout" )
    job.err = os.path.join( directory, "stderr" )

    commands = []
    if isinstance( job, JobGroup ):
      for task in range(1,job.tasks+1):
        command = "export SGE_TASK_ID=%d; %s" % (task, job.scriptPath)
        commands.append(command)
    else:
      commands.append(job.scriptPath)

    count = 0
    for command in commands:
      #print "# %s" % (command)
      #command += " 2>%s/%s.%d >%s/%s.%d" % (job.err, job.name, count, job.out, job.name, count)
      #subprocess.call(command)
      #print 'the command is', command
      errorCode = pool.apply_async(subprocess.call, (command,))
      errorCodes.append(errorCode)
      count += 1
    job.submitted = 1
  # wait for submitted jobs to finish
  pool.close()
  pool.join()
  # make sure all the jobs were successful
  for index, code in enumerate(errorCodes):
    #print index, code
    try:
      code.successful()
    except AssertionError, e:
      raise AssertionError("Job Failed to run: %s  %s" % (jobs[index], commands[index]), e.args)

def check_for_qsub():
  for directory in os.environ['PATH'].split(':'):
    if os.path.isfile(os.path.join(directory, "qsub")):
      return True
  return False

def submit_jobs( directory, jobs, use_grid_engine=True ):
  waiting = list(jobs)

  while len(waiting) > 0:
    # extract submittable jobs
    submittable = extract_submittable_jobs( waiting )

    # run those jobs
    if use_grid_engine:
      submit_safe_jobs( directory, submittable )
    else:
      run_safe_jobs( directory, submittable )

    # remove those from the waiting list
    map( waiting.remove, submittable )


def build_submission( directory, jobs, use_grid_engine=True ):
  # check to see if qsub exists
  if use_grid_engine: # make sure that SGE is available
    use_grid_engine = check_for_qsub()

  # build all necessary directories
  build_directories( directory )

  # build job scripts
  build_job_scripts( directory, jobs, use_grid_engine )

  # submit the jobs
  submit_jobs( directory, jobs, use_grid_engine )


def main():
  ''' map() a function against data as given in an input pickle file, saving result to a pickle.'''
  import os, optparse

  usage = "%prog [options] inputPickle outputPickle \n" + main.__doc__
  parser = optparse.OptionParser(usage)
  parser.add_option('--mode', dest='mode', type='string',
            help="""set the mode of operation.  Should be one of ['map', 'reduce']""")
  opts, args = parser.parse_args()
  inputPickleName, outputPickleName = args
  with open(inputPickleName) as inputPickleFile:
    with open(outputPickleName, 'wb') as outputPickleFile:
      os.environ = pickle.load(inputPickleFile)
      try:
        initializer = pickle.load(inputPickleFile)
        if opts.mode == 'map':
          initargs = pickle.load(inputPickleFile)
        else:
          initargs = None
      except Exception as e:
        # save and return the exception, to be debugged upstream
        pickle.dump(type(e)(('Error Loading initializer from pickle file',) + e.args), outputPickleFile)
        return
      else:
        if initializer is not None:
          if initargs is not None:
            initializer(initargs)
          else:
            initializer()
      try:
        func = pickle.load(inputPickleFile)
      except AttributeError as e:
        pickle.dump(e, outputPickleFile)  # Error loading
        return
      else:
        # If I didn't want error checking on each step, I'd do this at the top level:
        ## result = reduce(func, (data if type(data) != Job else pickle.load(data.outputFile) for data in curData), initializer)
        moreData = True
        curData = []
        while moreData:
          try:
            data = pickle.load(inputPickleFile)
          except EOFError:
            moreData = False
          else:
            if data.__class__ == Job:  # load the results from the previous Job
              try:
                data = pickle.load(data.outputFile)
              except Exception as e:
                pickle.dump(type(e)(("Error loading previous Job's results! ",) + e.args), outputPickleFile)
                return
            curData.append(data)
            if opts.mode == 'map':
              try:
                result = func(curData.pop(0))
              except Exception as e:
                # save and return the exception, to be debugged upstream
                pickle.dump(type(e)(('Error in map function: ',) + e.args), outputPickleFile)
                return
              else:
                pickle.dump(result, outputPickleFile)
            else:
              # reduce curData
              result = reduce(func, curData[:2], initializer)
              initializer = None
              del curData[:2]
              curData.insert(0, result)


if __name__ == '__main__':
  main()
