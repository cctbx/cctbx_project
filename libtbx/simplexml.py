from __future__ import absolute_import, division, print_function
class SimpleNode:
  def __init__(self,tag,contents='',indent=True):
    self.tag=tag
    self.m_attributes=[]
    self.children=[]
    self.content=contents
    self.indent_f=indent
  def attribute(self,key,value):
    self.m_attributes.append((key,value))
    return self
  def child(self,node):
    self.children.append(node)
    return node
  def contents(self,c):
    self.content=c
    return self
  def emit(self,channel,indent=0):
    if self.indent_f==False: indent=0
    attrs = [' %s="%s"'%(item[0],item[1]) for item in self.m_attributes]
    all_attrs = "".join(attrs)

    if self.content!='' and len(self.content)<80:
      print("%s<%s%s>%s</%s>"%(' '*indent,self.tag,all_attrs,self.content,self.tag), file=channel)
      return
    print("%s<%s%s>"%(' '*indent,self.tag,all_attrs), end=' ', file=channel)
    if self.content!='' and len(self.content)>=80:
      print(file=channel)
      print(self.content, file=channel)
    else: print(file=channel)
    for item in self.children:
      item.emit(channel,indent=indent+2)
    if len(self.children)>0:
      print("%s</%s>"%(' '*indent,self.tag), file=channel)
    else:
      channel.seek(channel.tell()-3)
      print("/>", file=channel)
