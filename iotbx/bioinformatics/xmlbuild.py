"""
Tools for working with xml
"""

from __future__ import absolute_import, division, print_function

class Node(object):
  """
  A node that allows setting arbitrary attributes
  """

def dummy(value):

  return value

# Extractions
class Value(object):
  """
  Extract text from elem
  """

  def __init__(self, attribute, conversion = dummy):

    self.attribute = attribute
    self.conversion = conversion


  def __call__(self, elem, node):

    setattr( node, self.attribute, self.conversion( elem.text ) )


class Text(object):
  """
  Extract text
  """

  def __init__(self, attribute, tag, conversion = dummy):

    self.attribute = attribute
    self.tag = tag
    self.conversion = conversion


  def __call__(self, elem, node):

    child = elem.find( self.tag )

    if child is None:
      raise RuntimeError("Element '%s' has no child '%s'" % (elem.tag, self.tag ))

    setattr( node, self.attribute, self.conversion( child.text ) )


class TextWithDefault(object):
  """
  Extract text or fill in default if not present
  """

  def __init__(self, attribute, tag, default = None, conversion = dummy):

    self.attribute = attribute
    self.tag = tag
    self.default = default
    self.conversion = conversion


  def __call__(self, elem, node):

    child = elem.find( self.tag )

    if child is None:
      value = self.default

    else:
      value = self.conversion( child.text )

    setattr( node, self.attribute, value )


class Attribute(object):
  """
  Extract an attribute
  """

  def __init__(self, attribute, name, conversion = dummy):

    self.attribute = attribute
    self.name = name
    self.conversion = conversion


  def __call__(self, elem, node):

    res = elem.get( self.name )

    if res is None:
      raise RuntimeError("Element '%s' has no attribute '%s'" % (
        elem.tag,
        self.name,
        ))

    setattr( node, self.attribute, self.conversion( res ) )


class AttributeWithDefault(object):
  """
  Extract an attribute or fill in a default value
  """

  def __init__(self, attribute, name, default, conversion = dummy):

    self.attribute = attribute
    self.name = name
    self.default = default
    self.conversion = conversion


  def __call__(self, elem, node):

    res = elem.get( self.name, self.default )
    setattr( node, self.attribute, self.conversion( res ) )


class Group(object):
  """
  Adds extra node to group items with a logical connection
  """

  def __init__(self, attribute, extractions):

    self.attribute = attribute
    self.extractions = extractions


  def __call__(self, elem, node):

    child = Node()
    setattr( node, self.attribute, child )

    for processing in self.extractions:
      processing( elem = elem, node = child )


# Nodes
class Single(object):
  """
  Extract data from an ElementTree Element into a Node
  """

  def __init__(self, child_data_tagged = {}, extractions = []):

    self.child_data_tagged = child_data_tagged
    self.extractions = extractions


  def __call__(self, iterstream, endtag):

    unseen = set( self.child_data_tagged )
    node = Node()

    for ( event, elem ) in iterstream:
      if elem.tag == endtag:
        assert event == "end"

        for processing in self.extractions:
          processing( elem = elem, node = node )

        elem.clear()
        break

      elif event == "start" and elem.tag in self.child_data_tagged:
        unseen.remove( elem.tag )
        ( attribute, processor ) = self.child_data_tagged[ elem.tag ]
        child = processor( iterstream = iterstream, endtag = elem.tag )
        setattr( node, attribute, child )

    if unseen:
      raise RuntimeError("Missing children for element '%s': %s" % (
        endtag,
        ", ".join( "'%s'" % t for t in unseen ),
        ))

    return node


class Multiple(object):
  """
  Extract data from an ElementTree Element into a child Node
  """

  def __init__(self, tag, processor):

    self.tag = tag
    self.processor = processor


  def __call__(self, iterstream, endtag):

    array = []

    for ( event, elem ) in iterstream:
      if event == "end":
        assert elem.tag == endtag
        elem.clear()
        break

      elif event == "start":
        assert elem.tag == self.tag
        child = self.processor( iterstream = iterstream, endtag = self.tag )
        array.append( child )

    return array


class DataAttribute(object):
  """
  A node that is converted directly to data, i.e. no children
  """

  def __init__(self, name, conversion = dummy):

    self.name = name
    self.conversion = conversion


  def __call__(self, iterstream, endtag):

    for ( event, elem ) in iterstream:
      if event == "end" and elem.tag == endtag:
        raw = elem.get( self.name )

        if raw is None:
          raise RuntimeError("%s has no %s attribute" % ( elem.tag, self.name ))

        result = self.conversion( raw )
        elem.clear()
        return result


# Parser
class Parser(object):

  def __init__(self, tag, builder, restype, cElementTree = False):

    self.tag = tag
    self.builder = builder
    self.restype = restype

    if cElementTree:
      import xml.etree.cElementTree
      self.module = xml.etree.cElementTree

    else:
      import xml.etree.ElementTree
      self.module = xml.etree.ElementTree


  def __call__(self, source):

    iterstream = self.module.iterparse( source, events = ( "start", "end" ) )

    for ( event, elem ) in iterstream:
      if elem.tag == self.tag and event == "start":
        return self.restype(
          root = self.builder( iterstream = iterstream, endtag = self.tag )
          )

    else:
      raise RuntimeError("Start tag %s not found" % self.tag)

