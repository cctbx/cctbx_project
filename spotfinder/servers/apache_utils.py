from __future__ import absolute_import, division, print_function
from libtbx.simplexml import SimpleNode

class LongLineSimpleNode(SimpleNode):

  def __init__(self,*args,**kwargs):
    SimpleNode.__init__(self,*args,**kwargs)

  def emit(self,channel,indent=0):
    if self.indent_f==False: indent=0
    attrs = [' %s="%s"'%(item[0],item[1]) for item in self.m_attributes]
    all_attrs = "".join(attrs)

    if self.content!='':
      print("%s<%s%s>%s</%s>"%(' '*indent,self.tag,all_attrs,self.content,self.tag), file=channel)
      return
    print("%s<%s%s>"%(' '*indent,self.tag,all_attrs), end=' ', file=channel)
    print(file=channel)
    for item in self.children:
      item.emit(channel,indent=indent+2)
    if len(self.children)>0:
      print("%s</%s>"%(' '*indent,self.tag), file=channel)
    else:
      channel.seek(channel.tell()-2)
      print("/>", file=channel)
