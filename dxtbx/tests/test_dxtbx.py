from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.sweep_filenames import template_regex

def test_template_regex():
  questions_answers = {
      'foo_bar_001.img':'foo_bar_###.img',
      'foo_bar001.img':'foo_bar###.img',
      'foo_bar_1.8A_001.img':'foo_bar_1.8A_###.img',
      'foo_bar.001':'foo_bar.###',
      'foo_bar_001.img1000':'foo_bar_###.img1000',
      'foo_bar_00001.img':'foo_bar_#####.img'
      }

  for filename in questions_answers:
    answer = template_regex(filename)
    assert answer[0] == questions_answers[filename]
