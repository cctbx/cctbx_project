from __future__ import absolute_import, division, print_function

import os

import pytest

tests_directory = os.path.dirname(os.path.realpath(__file__))

@pytest.fixture
def testsdir():
  return tests_directory
