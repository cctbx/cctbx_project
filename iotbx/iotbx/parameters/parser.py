from iotbx import parameters

def collect_values(word_stack, last_word, stop_token):
  result = []
  while (len(word_stack) > 0):
    word = word_stack.pop()
    if ((last_word.value == "\\" and last_word.quote_token is None)
        or word.quote_token is not None):
      result.append(word)
    elif (word.line_number != last_word.line_number
          or (stop_token is not None
              and word.value == stop_token and word.quote_token is None)):
      word_stack.push_back(word)
      break
    elif (word.value != "\\" or word.quote_token is not None):
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
      if (not parameters.is_standard_identifier(word.value)):
        word.raise_syntax_error("improper table name ")
      table = parameters.table(name=word.value, row_names=[], row_objects=[])
      while True:
        word = word_stack.pop_unquoted()
        if (word.value[0] != "."): break
        if (not table.has_attribute_with_name(word.value[1:])):
          raise RuntimeError(
            "Unexpected table attribute: %s (input line %d)" % (
              word.value, word.line_number))
        table.assign_attribute(
          name=word.value[1:],
          value_words=collect_values(word_stack, word, stop_token))
      word.assert_expected("{")
      word = word_stack.pop_unquoted()
      while (word.value != "}"):
        row_name = None
        if (word.value != "{"):
          row_name = word.value
          if (not parameters.is_standard_identifier(row_name)):
            word.raise_syntax_error("improper table row name ")
          word = word_stack.pop_unquoted()
          word.assert_expected("{")
        table.add_row(
          name=row_name,
          objects=collect_objects(
            word_stack=word_stack,
            definition_type_names=definition_type_names,
            stop_token="}"))
        word = word_stack.pop_unquoted()
      objects.append(table)
    else:
      lead_word = word
      if (lead_word.value == "{"):
        raise RuntimeError(
          'Syntax error: unexpected "{" (input line %d)' % (
            lead_word.line_number))
      word = word_stack.pop()
      if (word.quote_token is None and word.value == "{"):
        if (not parameters.is_standard_identifier(lead_word.value)):
          lead_word.raise_syntax_error("improper scope name ")
        objects.append(parameters.scope(
          name=lead_word.value,
          objects=collect_objects(
            word_stack=word_stack,
            definition_type_names=definition_type_names,
            stop_token="}")))
      else:
        word_stack.push_back(word)
        if (lead_word.value[0] != "."):
          if (not parameters.is_standard_identifier(lead_word.value)):
            lead_word.raise_syntax_error("improper definition name ")
          objects.append(parameters.definition(
            name=lead_word.value,
            values=collect_values(word_stack, lead_word, stop_token)))
        else:
          if (len(objects) == 0
              or not isinstance(objects[-1], parameters.definition)
              or not objects[-1].has_attribute_with_name(lead_word.value[1:])):
            raise RuntimeError("Unexpected attribute: %s (input line %d)" % (
              lead_word.value, lead_word.line_number))
          objects[-1].assign_attribute(
            name=lead_word.value[1:],
            value_words=collect_values(word_stack, lead_word, stop_token),
            type_names=definition_type_names)
  if (stop_token is not None):
    raise RuntimeError('Syntax error: missing "%s".' % stop_token)
  return objects
