import generate_reductions
import generate_std_imports
import generate_operator_functors
import generate_operator_traits_builtin
import generate_algebras
import generate_apply

def run():
  generate_reductions.run()
  generate_std_imports.run()
  generate_operator_functors.run()
  generate_operator_traits_builtin.run()
  generate_algebras.run()
  generate_apply.run()

if (__name__ == "__main__"):
  run()
