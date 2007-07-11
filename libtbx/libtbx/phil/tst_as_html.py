from libtbx import phil
import libtbx.phil.as_html

test_params = phil.parse(
  input_string="""\
par00 = 0
  .help = abc00
  .type = float
par01 = False
  .help = abc01
  .type = bool
scope1
  .help = "scope1 help"
{
  par10 = 1
    .help = abc11
    .type = float
  par11 = True
    .type = bool
  scope20
    .help = "very long scope20 help very long scope20 help very long scope20"
            "help very long scope20 help very long scope20 help very long"
            "scope20 help"
  {
    par20 = 2
      .type = float
    par21 = False
      .type = bool
    scope3 {
      par30 = 3
        .type = float
      par31 = True
        .type = bool
    }
  }
  scope21 {
    par21 = 2
      .help = "very long par21 help very long par21 help very long par21 help"
              "very long par21 help very long par21 help very long par21 help"
              "very long par21 help very long par21 help very long par21 help"
      .type = float
    par21 = False
      .type = bool
    scope3 {
      par30 = 3
        .type = float
      par31 = True
        .type = bool
    }
  }
  par12 = 10
    .type = float
}
par02 = 0
  .type = float
par03 = False
  .type = bool
""",
  process_includes=True)

def run():
  phil.as_html.run(phil_object=test_params, log=open("tst_as_html.html", "w"))
  print "OK"

if (__name__ == "__main__" ):
  run()
