---------------------------------------------------------------------
libtbx.easy_mp: functions for painless parallelization of Python code
---------------------------------------------------------------------


pool_map: parallel map() function for Unix/Linux systems
========================================================

``pool_map`` is essentially a simple replacement for the built-in map function, and is fast enough to be suitable for parallelization across arbitrarily large
shared-memory systems even for relatively short-running methods.  Although it
is based on the ``map`` method of the ``multiprocessing.Pool`` class, the
underlying implementation has been modified to take advantage of the internal
behavior.  On systems that support the ``os.fork()`` call, it can embed a
reference to a "fixed" function that is preserved in the child processes,
bypassing all ``pickle`` calls.  This means that unpickleable objects can be
embedded in a callable which is given as the ``fixed_func`` argument.
``easy_mp`` is therefore suitable for parallelizing existing code that would
otherwise require extensive refactoring (for instance, the weight optimization
in phenix.refine).

Note that since ``pool_map`` depends on ``os.fork()``, it does not support
parallelization on Windows systems.  However, it will simply run in serial on
Windows, allowing it to be used in a platform-independent manner.

.. autofunction:: libtbx.easy_mp.pool_map
.. autofunction:: libtbx.easy_mp.detect_problem
.. autofunction:: libtbx.easy_mp.enable_multiprocessing_if_possible
.. autofunction:: libtbx.easy_mp.get_processes


parallel_map: parallel map() function for multiple architectures
================================================================

``parallel_map`` is a more general replacement for ``map()``, and supports
Windows systems in addition to Unix/Linux.  It also enables use of managed
clusters via several common queuing systems such as Sun Grid Engine or PBS.
Because of the higher overhead involved, it is more suitable for functions
that take minutes or longer to run, but the use of cluster resources means that
it can handle much larger operations efficiently.

.. autofunction:: libtbx.easy_mp.parallel_map
