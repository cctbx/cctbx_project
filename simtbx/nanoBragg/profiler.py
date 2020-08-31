#!/usr/bin/env python
# -*- coding: utf-8 -*-


import importlib
TIMEMORY_LOADER    = importlib.find_loader("timemory")
TIMEMORY_AVAILABLE = TIMEMORY_LOADER is not None


#
# Singleton class -- to store profiler
#

class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):

        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)

        return cls._instances[cls]



class Profiler(object, metaclass=Singleton):

    def __init__(self):
        self._enabled = False


    @property
    def enabled(self):
        return self._enabled


    @enabled.setter
    def enabled(self, val):
        self._enabled = val


    @property
    def profile_decorator(self):

        def identity2(*args, **kwargs):
            def identity(f):
                return f
            return identity
 

        if self.enabled:
            from timemory.profiler import profile
            return profile
        else:
            return identity2



def enable():
    if TIMEMORY_AVAILABLE:
        Profiler().enabled = True
    else:
        raise RuntimeError("Cannot enable profiling => Timemory not available")



def disable():
    Profiler().enables = False
/
