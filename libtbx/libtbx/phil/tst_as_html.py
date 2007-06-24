import iotbx.phil
from libtbx.phil import as_html

test_params = iotbx.phil.parse("""\
par00 = 0
  .type=float
  .help = abc00
par01 = False
  .type=bool
  .help = abc01
scope1 {
  par10 = 1
    .type=float
    .help = abc11
  par11 = True
    .type=bool
  scope20 {
    par20 = 2
      .type=float
    par21 = False
      .type=bool
    scope3 {
      par30 = 3
        .type=float
      par31 = True
        .type=bool
    }
  }
  scope21 {
    par21 = 2
      .type=float
    par21 = False
      .type=bool
    scope3 {
      par30 = 3
        .type=float
      par31 = True
        .type=bool
    }
  }
par12=10
  .type=float
}
par02 = 0
  .type=float
par03 = False
  .type=bool
""", process_includes=True)


def run():
  log = open("example.html", "w")
  try:
    #import phenix.command_line.lsq_superpose_pdbs
    from  phenix.command_line import lsq_superpose_pdbs
    master_params=lsq_superpose_pdbs.master_params

  except:
    raise RuntimeError("Cannot import phenix.refinement.")
  as_html.run(phil_object = master_params,
              log         = log)

if (__name__ == "__main__" ):
  run()
