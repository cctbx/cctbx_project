from __future__ import absolute_import, division, print_function

import unittest
import time
from six.moves import range

def normal_function(arg, wait):

  time.sleep( wait )
  return 3 * arg


def raise_runtime_error(message, wait):

  time.sleep( wait )
  raise RuntimeError(message)


def raise_sorry(message, wait):

  time.sleep( wait )
  from libtbx.utils import Sorry
  raise Sorry(message)


def crash_python(wait):

  time.sleep( wait )
  import boost_adaptbx.boost.python as bp
  bp.ext.dereference_char_pointer(None)


class TestManager(unittest.TestCase):

  def run_normal_function(self, manager, value):

    self.assertTrue( manager.is_empty() )
    jobid = manager.submit( target = normal_function, args = ( value, 0.05 ) )
    self.assertFalse( manager.is_empty() )
    result = next(manager.results())
    self.assertEqual( result[0], jobid )
    self.assertEqual( normal_function( arg = value, wait = 0 ), result[1]() )
    self.assertTrue( manager.is_empty() )


  def run_raise_runtime_error(self, manager, message):

    self.assertTrue( manager.is_empty() )
    jobid = manager.submit( target = raise_runtime_error, args = ( message, 0.05 ) )
    self.assertFalse( manager.is_empty() )
    result = next(manager.results())
    self.assertEqual( result[0], jobid )
    self.assertRaises( RuntimeError, result[1] )
    self.assertTrue( manager.is_empty() )


  def run_raise_sorry(self, manager, message):

    self.assertTrue( manager.is_empty() )
    jobid = manager.submit( target = raise_sorry, args = ( message, 0.05 ) )
    self.assertFalse( manager.is_empty() )
    result = next(manager.results())
    self.assertEqual( result[0], jobid )

    from libtbx.utils import Sorry
    self.assertRaises( Sorry, result[1] )
    self.assertTrue( manager.is_empty() )


  def run_crash_python(self, manager):

    self.assertTrue( manager.is_empty() )
    jobid = manager.submit( target = crash_python, args = ( 0.05, ) )
    self.assertFalse( manager.is_empty() )
    result = next(manager.results())
    self.assertEqual( result[0], jobid )

    self.assertRaises( RuntimeError, result[1] )
    self.assertTrue( manager.is_empty() )


  def run_tests(self, creator, perform_crash_test = False):

    from libtbx.scheduling import holder

    with holder( creator = creator ) as manager:
      for i in range( 5 ):
        self.run_normal_function( manager = manager, value = i )

      self.run_raise_runtime_error( manager = manager, message = "runtime_error message" )
      self.run_raise_sorry( manager = manager, message = "sorry message" )

      if perform_crash_test:
        self.run_crash_python( manager = manager )

      for i in range( 5, 10 ):
        self.run_normal_function( manager = manager, value = i )

      manager.join()


  def test_mainthread(self):

    from libtbx.scheduling import mainthread
    self.run_tests( creator = mainthread.creator )


  def test_job_scheduler_threading(self):

    from libtbx.scheduling import job_scheduler
    from libtbx.scheduling import thread_handler
    import threading

    creator = job_scheduler.creator(
      job_factory = threading.Thread,
      queue_factory = thread_handler.qfactory,
      capacity = job_scheduler.limited( njobs = 1 ),
      )
    self.run_tests( creator = creator )


  def test_job_scheduler_multiprocessing(self):

    from libtbx.scheduling import job_scheduler
    from libtbx.scheduling import mp_handler

    creator = job_scheduler.creator(
      job_factory = mp_handler.stderr_capturing_process,
      queue_factory = mp_handler.fifo_qfactory,
      capacity = job_scheduler.limited( njobs = 1 ),
      )
    import sys
    self.run_tests(
      creator = creator,
      perform_crash_test = ( sys.platform != "win32" ),
      )


  def test_process_pool_threading(self):

    from libtbx.scheduling import process_pool
    from libtbx.scheduling import thread_handler
    import threading

    creator = process_pool.creator(
      job_factory = threading.Thread,
      inq_factory = thread_handler.qfactory,
      outq_factory = thread_handler.qfactory,
      autoscaling = process_pool.constant_capacity( capacity = 1 ),
      lifecycle = process_pool.unlimited,
      )
    self.run_tests( creator = creator )


  def test_process_pool_multiprocessing(self):

    from libtbx.scheduling import process_pool
    from libtbx.scheduling import mp_handler

    creator = process_pool.creator(
      job_factory = mp_handler.stderr_capturing_process,
      inq_factory = mp_handler.fifo_qfactory,
      outq_factory = mp_handler.fifo_qfactory,
      autoscaling = process_pool.constant_capacity( capacity = 1 ),
      lifecycle = process_pool.unlimited,
      )
    import sys
    self.run_tests(
      creator = creator,
      perform_crash_test = ( sys.platform != "win32" ),
      )


suite_manager = unittest.TestLoader().loadTestsFromTestCase(
  TestManager
  )


alltests = unittest.TestSuite(
  [
    suite_manager,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )
