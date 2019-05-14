from __future__ import absolute_import, division, print_function
import scitbx.lbfgs
from scitbx.array_family import flex

class refinery:

  def __init__(self):
    print("refinery")
    self.x = flex.double([0])
    scitbx.lbfgs.run(target_evaluator=self)

  def compute_functional_and_gradients(self):
    print("compute_functional_and_gradients")
    f = 0
    g = flex.double([0])
    return f, g

  def callback_after_step(self, minimizer):
    print("callback_after_step")

def run():
  refinery()
  print("OK")

if (__name__ == "__main__"):
  run()
