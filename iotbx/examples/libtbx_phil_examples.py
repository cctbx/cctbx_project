"""Examples of use of the PHIL command syntax"""
from __future__ import absolute_import, division, print_function
from six.moves import range

if (__name__ == "__main__"):

  # ---- line 197 -------------------------------------------------------------

  from libtbx.phil import parse

  master_phil = parse("""
    minimization.input {
      file_name = None
        .type = path
      label = None
        .type = str
    }
    """)

  user_phil = parse("""
    minimization.input {
      file_name = experiment.dat
    }
    """)

  command_line_phil = parse(
    "minimization.input.label=set2")

  working_phil = master_phil.fetch(
    sources=[user_phil, command_line_phil])
  working_phil.show()

  # ---- line 241 -------------------------------------------------------------

  argument_interpreter = master_phil.command_line_argument_interpreter(
    home_scope="minimization")

  command_line_phil = argument_interpreter.process(
    arg="minimization.input.label=set2")

  # ---- line 281 -------------------------------------------------------------

  working_params = working_phil.extract()

  # ---- line 289 -------------------------------------------------------------

  print(working_params.minimization.input.file_name)
  print(working_params.minimization.input.label)

  # ---- line 311 -------------------------------------------------------------

  working_params.minimization.input.label = "set3"
  modified_phil = master_phil.format(python_object=working_params)
  modified_phil.show()

  # ---- line 343 -------------------------------------------------------------

  master_phil = parse("""
    minimization.input {
      file_name = None
        .type = path
        .multiple = True
    }
    """)

  # ---- line 357 -------------------------------------------------------------

  user_phil = parse("""
    minimization.input {
      file_name = experiment1.dat
      file_name = experiment2.dat
      file_name = experiment3.dat
    }
    """)

  # ---- line 377 -------------------------------------------------------------

  working_params = master_phil.fetch(source=user_phil).extract()
  print(working_params.minimization.input.file_name)

  # ---- line 390 -------------------------------------------------------------

  master_phil = parse("""
    minimization {
      input
        .multiple = True
      {
        file_name = None
          .type = path
        label = None
          .type = str
      }
    }
    """)

  # ---- line 409 -------------------------------------------------------------

  user_phil = parse("""
    minimization {
      input {
        file_name = experiment1.dat
        label = set2
      }
      input {
        file_name = experiment2.dat
        label = set1
      }
    }
    """)

  # ---- line 428 -------------------------------------------------------------

  working_params = master_phil.fetch(source=user_phil).extract()
  for input in working_params.minimization.input:
    print(input.file_name)
    print(input.label)

  # ---- line 448 -------------------------------------------------------------

  master_phil = parse("""
    minimization {
      input
        .multiple = True
      {
        file_name = None
          .type = path
        label = None
          .type = str
          .multiple = True
      }
    }
    """)

  # ---- line 468 -------------------------------------------------------------

  user_phil = parse("""
    minimization {
      input {
        file_name = experiment1.dat
        label = set1
        label = set2
        label = set3
      }
      input {
        file_name = experiment2.dat
        label = set2
        label = set3
      }
    }
    """)

  # ---- line 490 -------------------------------------------------------------

  working_params = master_phil.fetch(source=user_phil).extract()
  for input in working_params.minimization.input:
    print(input.file_name)
    print(input.label)

  # ---- line 519 -------------------------------------------------------------

  master_phil = parse("""
    minimization.parameters {
      method = *bfgs conjugate_gradient
        .type = choice
      max_iterations = 10
        .type = int
    }
    """)

  user_phil = parse("""
    minimization.parameters {
      method = bfgs *conjugate_gradient
    }
    """)

  working_phil = master_phil.fetch(source=user_phil)
  diff_phil = master_phil.fetch_diff(source=working_phil)
  diff_phil.show()

  # ---- line 599 -------------------------------------------------------------

  var_phil = parse("""
    root_name = peak
    file_name = $root_name.mtz
    full_path = $HOME/$file_name
    related_file_name = $(root_name)_data.mtz
    message = "Reading $file_name"
    as_is = ' $file_name '
    """)
  var_phil.fetch(source=var_phil).show()

  # ---- line 649 -------------------------------------------------------------

  import libtbx.phil
  from libtbx.phil import tokenizer

  class upper_converters:

    phil_type = "upper"

    def __str__(self): return self.phil_type

    def from_words(self, words, master):
      s = libtbx.phil.str_from_words(words=words)
      if (s is None): return None
      return s.upper()

    def as_words(self, python_object, master):
      if (python_object is None):
        return [tokenizer.word(value="None")]
      return [tokenizer.word(value=python_object.upper())]

  converter_registry = libtbx.phil.extended_converter_registry(
    additional_converters=[upper_converters])

  # ---- line 678 -------------------------------------------------------------

  master_phil = parse("""
    value = None
      .type = upper
    """,
      converter_registry=converter_registry)
  user_phil = parse("value = extracted")
  working_params = master_phil.fetch(source=user_phil).extract()
  print(working_params.value)

  # ---- line 694 -------------------------------------------------------------

  working_params.value = "formatted"
  working_phil = master_phil.format(python_object=working_params)
  working_phil.show()

  # ---- line 727 -------------------------------------------------------------

  master_phil = parse("""
    random_integers = None
      .type = ints
    euler_angles = None
      .type = floats(size=3)
    unit_cell_parameters = None
      .type = floats(size_min=1, size_max=6)
    rotation_part = None
      .type = ints(size=9, value_min=-1, value_max=1)
    """)

  user_phil = parse("""
    random_integers = 3 18 5
    euler_angles = 10 -20 30
    unit_cell_parameters = 10,20,30
    rotation_part = "1,0,0;0,-1,0;0,0,-1"
    """)

  working_phil = master_phil.fetch(source=user_phil)
  working_phil.show()
  print()
  working_params = working_phil.extract()
  print(working_params.random_integers)
  print(working_params.euler_angles)
  print(working_params.unit_cell_parameters)
  print(working_params.rotation_part)
  print()
  working_phil = master_phil.format(python_object=working_params)
  working_phil.show()

  # ---- line 805 -------------------------------------------------------------

  master_phil = parse("""
    gender = male female
      .type = choice
    favorite_sweets = ice_cream chocolate candy_cane cookies
      .type = choice(multi=True)
    """)

  jims_choices = parse("""
    gender = *male female
    favorite_sweets = *ice_cream chocolate candy_cane *cookies
    """)

  jims_phil = master_phil.fetch(source=jims_choices)
  jims_phil.show()
  jims_params = jims_phil.extract()
  print(jims_params.gender, jims_params.favorite_sweets)

  # ---- line 845 -------------------------------------------------------------

  ignorant_choices = parse("""
    gender = male female
    favorite_sweets = ice_cream chocolate candy_cane cookies
    """)

  ignorant_params = master_phil.fetch(source=ignorant_choices).extract()
  print(ignorant_params.gender, ignorant_params.favorite_sweets)

  # ---- line 880 -------------------------------------------------------------

  greedy_choices = parse("""
    favorite_sweets=ice_cream+chocolate+cookies
    """)

  greedy_params = master_phil.fetch(source=greedy_choices).extract()
  print(greedy_params.favorite_sweets)

  # ---- line 898 -------------------------------------------------------------

  no_thanks_choices = parse("""
    favorite_sweets=None
    """)

  no_thanks_params = master_phil.fetch(source=no_thanks_choices).extract()
  print(no_thanks_params.favorite_sweets)

  # ---- line 920 -------------------------------------------------------------

  master_phil = parse("""
    minimization.input {
      file_name = None
        .type = path
    }
    minimization.parameters {
      max_iterations = 10
        .type = int
    }
    """)

  user_phil = parse("""
    minimization.input.file_name = experiment.dat
    minimization.parameters.max_iterations = 5
    """)

  working_params = master_phil.fetch(source=user_phil).extract()
  print(working_params)
  print(working_params.minimization.input.file_name)
  print(working_params.minimization.parameters.max_iterations)

  # ---- line 956 -------------------------------------------------------------

  print(working_params.minimization.input.__phil_path__())
  print(working_params.minimization.parameters.__phil_path__())

  # ---- line 1014 ------------------------------------------------------------

  master_phil = parse("""
    plot
      .multiple = True
    {
      style = line bar pie_chart
        .type=choice
      title = None
        .type = str
    }
    plot {
      style = line
      title = Line plot (default in master)
    }
    """)

  user_phil = parse("""
    plot {
      style = bar
      title = Bar plot (provided by user)
    }
    """)

  working_phil = master_phil.fetch(source=user_phil)
  working_phil.show()

  # ---- line 1056 ------------------------------------------------------------

  working_params = working_phil.extract()
  print(working_params.plot)

  # ---- line 1073 ------------------------------------------------------------

  master_phil = parse("""
    plot
      .multiple = True
      .optional = False
    {
      style = line bar pie_chart
        .type=choice
      title = None
        .type = str
    }
    plot {
      style = line
      title = Line plot (default in master)
    }
    """)

  # ---- line 1096 ------------------------------------------------------------

  working_phil = master_phil.fetch(source=user_phil)
  working_phil.show()
  print(working_phil.extract().plot)

  # ---- line 1144 ------------------------------------------------------------

  master_phil = parse("""
    input {
      file_name = None
        .type = path
    }
    """)

  user_phil = parse("""
    input {
      file_name = experiment.dat
      label = set1
      lable = set2
    }
    """)

  working_phil, unused = master_phil.fetch(
    source=user_phil, track_unused_definitions=True)
  working_phil.show()
  for object_locator in unused:
    print("unused:", object_locator)

  # ---- line 1189 ------------------------------------------------------------

  phil_scope = parse("""
     quick .multiple=true;.optional=false{and=very;.type=str;dirty=use only on command-lines, please!;.type=str}
     """)

  phil_scope.show(attributes_level=2)

  # ---- line 1232 ------------------------------------------------------------

  master_phil = parse("""
    !input {
      file_name = None
        .type = path
        .multiple = True
    }
    """)
  master_phil.show()

  # ---- line 1254 ------------------------------------------------------------

  user_phil = parse("""
    input.file_name = experiment.dat
    """)
  print(len(master_phil.fetch(source=user_phil).as_str()))

  # ---- line 1316 ------------------------------------------------------------

  master_phil = parse("""
    minimization {
      input
        .help = "File names and data labels."
        .multiple = True
      {
        file_name = None
          .type = path
        label = None
          .help = "A unique substring of the data label is sufficient."
          .type = str
      }
    }
    """)

  for attributes_level in range(4):
    master_phil.show(attributes_level=attributes_level)
