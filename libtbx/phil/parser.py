from __future__ import absolute_import, division, print_function
import libtbx.phil

def collect_assigned_words(word_iterator, lead_word):
  have_comment = False
  last_word = lead_word
  result = []
  while True:
    word = word_iterator.try_pop(settings_index=1)
    if (word is None): break
    if (not have_comment
        and word.quote_token is None
        and word.value in ["{", "}", ";", "#"]):
      if (word.value == ";"):
        break
      if (word.value != "#"):
        word_iterator.backup()
        break
      have_comment = True
    elif (word.quote_token is not None
        or (last_word.quote_token is None and last_word.value == "\\")):
      if (not have_comment): result.append(word)
    elif (word.line_number != last_word.line_number):
      word_iterator.backup()
      break
    elif (word.value != "\\" or word.quote_token is not None):
      if (not have_comment): result.append(word)
    last_word = word
  if (len(result) == 0):
    raise RuntimeError("Missing value for %s%s" % (
      str(lead_word), lead_word.where_str()))
  return result

def collect_objects(
      word_iterator,
      converter_registry,
      primary_id_generator,
      primary_parent_scope,
      stop_token=None,
      start_word=None,
      definition_converter_cache=None,
      scope_extract_call_proxy_cache=None):
  if (definition_converter_cache is None):
    definition_converter_cache = {}
  if (scope_extract_call_proxy_cache is None):
    scope_extract_call_proxy_cache = {}
  prev_line_number = 0
  active_definition = None
  while True:
    lead_word = word_iterator.try_pop_unquoted()
    if (lead_word is None): break
    if (lead_word.value == "#phil"
        and lead_word.line_number != prev_line_number):
      word = word_iterator.pop_unquoted()
      if (word.value == "__END__"): break
      if (word.value == "__ON__"): continue
      if (word.value != "__OFF__"):
        raise RuntimeError("Unknown: %s %s%s" % (
          lead_word.value, word.value, word.where_str()))
      if (word_iterator.scan_for_start(
            intro=lead_word.value, followups=["__END__", "__ON__"]) == 0):
        break
      continue
    if (stop_token is not None and lead_word.value == stop_token):
      return
    if (lead_word.value == "{"):
      raise RuntimeError(
        'Syntax error: unexpected "{"%s' % lead_word.where_str())
    if (lead_word.value[:1] == "!"):
      lead_word.value = lead_word.value[1:]
      is_disabled = True
    else:
      is_disabled = False
    prev_line_number = lead_word.line_number
    word = word_iterator.pop()
    if (word.quote_token is None
        and (   word.value == "{"
             or word.value[:1] == "."
             or word.value[:2] == "!.")):
      if (not libtbx.phil.is_standard_identifier(lead_word.value)):
        if (lead_word.value == ";"):
          lead_word.raise_syntax_error("unexpected ")
        else:
          lead_word.raise_syntax_error("improper scope name ")
      active_definition = None
      scope = libtbx.phil.scope(
        name=lead_word.value,
        primary_id=next(primary_id_generator),
        is_disabled=is_disabled,
        where_str=lead_word.where_str())
      while True:
        if (word.value == "{"):
          break
        if (word.value[:1] == "!"):
          word.value = word.value[1:]
          is_disabled = True
        else:
          is_disabled = False
        if (word.value[:1] != "."
            or not scope.has_attribute_with_name(word.value[1:])):
          raise RuntimeError(
            "Unexpected scope attribute: %s%s" % (
              word.value, word.where_str()))
        word_iterator.pop_unquoted().assert_expected("=")
        assigned_words = collect_assigned_words(word_iterator, word)
        if (not is_disabled):
          scope.assign_attribute(
            name=word.value[1:],
            words=assigned_words,
            scope_extract_call_proxy_cache=scope_extract_call_proxy_cache)
        word = word_iterator.pop_unquoted()
      collect_objects(
        word_iterator=word_iterator,
        converter_registry=converter_registry,
        primary_id_generator=primary_id_generator,
        primary_parent_scope=scope,
        stop_token="}",
        start_word=word,
        definition_converter_cache=definition_converter_cache,
        scope_extract_call_proxy_cache=scope_extract_call_proxy_cache)
      primary_parent_scope.adopt(scope)
    else:
      word_iterator.backup()
      if (lead_word.value[:1] != "."):
        if (not libtbx.phil.is_standard_identifier(lead_word.value)):
          if (lead_word.value == ";"):
            lead_word.raise_syntax_error("unexpected ")
          else:
            lead_word.raise_syntax_error("improper definition name ")
        if (lead_word.value != "include"):
          word_iterator.pop_unquoted().assert_expected("=")
        active_definition = libtbx.phil.definition(
          name=lead_word.value,
          words=collect_assigned_words(word_iterator, lead_word),
          primary_id=next(primary_id_generator),
          is_disabled=is_disabled,
          where_str=lead_word.where_str())
        primary_parent_scope.adopt(active_definition)
      else:
        if (active_definition is None
            or not active_definition.has_attribute_with_name(
                     lead_word.value[1:])):
          raise RuntimeError("Unexpected definition attribute: %s%s" % (
            lead_word.value, lead_word.where_str()))
        word = word_iterator.pop_unquoted()
        word.assert_expected("=")
        assigned_words = collect_assigned_words(word_iterator, lead_word)
        if (not is_disabled):
          active_definition.assign_attribute(
            name=lead_word.value[1:],
            words=assigned_words,
            converter_registry=converter_registry,
            converter_cache=definition_converter_cache)
  if (stop_token is not None):
    if (start_word is None):
      where = ""
    else:
      where = start_word.where()
    if (where == ""):
      raise RuntimeError('Syntax error: missing "%s"' % stop_token)
    raise RuntimeError('Syntax error: no matching "%s" for "%s" at %s' %
      (stop_token, str(start_word), where))
