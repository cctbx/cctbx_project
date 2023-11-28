from .molstar import MolstarViewer

def is_available():
  module_dir = __file__ 
  print(f"Looking for Molstar in {module_dir}")
  return True