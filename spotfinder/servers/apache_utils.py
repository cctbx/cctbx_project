from libtbx.simplexml import SimpleNode

class LongLineSimpleNode(SimpleNode):

  def __init__(self,*args,**kwargs):
    SimpleNode.__init__(self,*args,**kwargs)

  def emit(self,channel,indent=0):
    if self.indent_f==False: indent=0
    attrs = [' %s="%s"'%(item[0],item[1]) for item in self.m_attributes]
    all_attrs = "".join(attrs)

    if self.content!='':
      print >>channel,"%s<%s%s>%s</%s>"%(' '*indent,self.tag,all_attrs,self.content,self.tag)
      return
    print >>channel,"%s<%s%s>"%(' '*indent,self.tag,all_attrs),
    print >>channel
    for item in self.children:
      item.emit(channel,indent=indent+2)
    if len(self.children)>0:
      print >>channel,"%s</%s>"%(' '*indent,self.tag)
    else:
      channel.seek(channel.tell()-2)
      print >>channel,"/>"
