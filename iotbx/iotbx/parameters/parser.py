from iotbx import parameters
from iotbx import simple_tokenizer
from cStringIO import StringIO

def collect_values(word_stack, last_word):
  result = []
  while (len(word_stack) > 0):
    word = word_stack.pop()
    if (last_word.value == "\\"):
      result.append(word)
    elif (word.line_number != last_word.line_number):
      word_stack.push_back(word)
      break
    elif (word.quote_char is not None or word.value != "\\"):
      result.append(word)
    last_word = word
  return result

def collect_objects(word_stack, definition_type_names, stop_token=None):
  objects = []
  while (len(word_stack) > 0):
    word = word_stack.pop_unquoted()
    if (stop_token is not None and word.value == stop_token):
      return objects
    if (word.value == "table"):
      word = word_stack.pop_unquoted()
      table = parameters.table(name=word.value)
      while True:
        word = word_stack.pop_unquoted()
        if (word.value[0] != "."): break
        if (not table.has_attribute_with_name(word.value[1:])):
          raise RuntimeError(
            "Unexpected table attribute: %s (input line %d)" % (
              word.value, word.line_number))
        table.assign_attribute(
          name=word.value[1:],
          value_words=collect_values(word_stack, word))
      assert word.value == "{"
      word = word_stack.pop_unquoted()
      while True:
        row_name = None
        if (word.value != "{"):
          row_name = word.value
          word = word_stack.pop_unquoted()
          assert word.value == "{"
        table.add_row(
          name=row_name,
          objects=collect_objects(
            word_stack=word_stack,
            definition_type_names=definition_type_names,
            stop_token="}"))
        word = word_stack.pop_unquoted()
        if (word.value == "}"): break
      objects.append(table)
    else:
      lead_word = word
      word = word_stack.pop()
      if (word.quote_char is None and word.value == "{"):
        objects.append(parameters.scope(
          name=lead_word.value,
          objects=collect_objects(
            word_stack=word_stack,
            definition_type_names=definition_type_names,
            stop_token="}")))
      else:
        word_stack.push_back(word)
        if (lead_word.value[0] != "."):
          objects.append(parameters.definition(
            name=lead_word.value,
            values=collect_values(word_stack, lead_word)))
        else:
          if (len(objects) == 0
              or not isinstance(objects[-1], parameters.definition)
              or not objects[-1].has_attribute_with_name(lead_word.value[1:])):
            raise RuntimeError("Unexpected attribute: %s (input line %d)" % (
              lead_word.value, lead_word.line_number))
          objects[-1].assign_attribute(
            name=lead_word.value[1:],
            value_words=collect_values(word_stack, lead_word),
            type_names=definition_type_names)
  if (stop_token is not None):
    raise RuntimeError('No matching "%s".' % stop_token)
  return objects

def run():
  import sys
  input_string = open(sys.argv[1]).read()
  word_stack = simple_tokenizer.as_word_stack(
    input_string=input_string,
    contiguous_word_characters="")
  objects = collect_objects(
    word_stack=word_stack,
    definition_type_names=parameters.default_definition_type_names)
  if (1):
    out=sys.stdout
  else:
    out=StringIO()
  for object in objects:
    object.show(out=out, prefix="", attributes_level=1)

run()
