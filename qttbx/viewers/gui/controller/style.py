from .controller import Controller
from ..state.style import Style
from ..state.ref import SelectionRef, MapRef, ModelRef, RestraintsRef

class ModelStyleController(Controller):
  """
  Manage styles. Tightly coupled to viewer.
  Expects viewer to be parent
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.state.signals.style_change.connect(self.apply_from_json)

  def transition_color(self,ref,style):
    self.parent.color_ref(ref,style.color)

  def transition_visible(self,ref, new_style):
    print("transition_visible() called...")

    if new_style.visible == False:
      self.parent.hide_ref(ref)
    else:
      self.parent.show_ref(ref)


    # elif isinstance(ref,SelectionRef):

    # print("New style:")
    # print(new_style.to_json(indent=2))

    # if new_style.visible == False:
    #   #self.parent.set_representation(ref.id,[])
    #   self.parent.set_visibility(ref.id,False,queue=True)
    # else:
    #   # last_reps = ref.last_visible_style.representations
    #   # self.parent.set_representation(ref.id,last_reps)
    #   base_rep = set(new_style.representation).intersection(set(oldstyle.representation))
    #   if len(base_rep)>0:
    #     self.parent.set_visibility(ref.id,True,queue=True)
    #   for rep_name in set(new_style.representation):#-set(base_rep):
    #     self.parent.set_representation(ref,rep_name,queue=True)

    #   if new_style.color is not None:
    #     self.transition_color(ref,new_style,queue=True)
    #   self.parent.reset_camera(queue=True)
    #   #self.parent.clear_selection(queue=True)

    # if not queue:
    #   self.parent.execute_command_queue()
      

  
  def transition_representation(self,ref,style,queue=False):
    new_style = style
    old_style = ref.style
    for rep_name in set(new_style.representation + old_style.representation):
      if rep_name in new_style.representation and rep_name not in old_style.representation:
        self.parent.show_representation(ref,rep_name)
      elif rep_name in old_style.representation and rep_name not in new_style.representation:
        self.parent.hide_representation(ref,rep_name)
    
  # Generic apply of a style to a ref. Uses above transition_ functions
  def apply_from_json(self,json_str):
    style = Style.from_json(json_str)
    self.apply(style)

  def apply(self, style):
    ref = self.state.references[style.ref_id]
    if isinstance(ref,(SelectionRef,RestraintsRef,ModelRef)): #TODO: Replace with 'model-like'
      prev_dict = ref.style.to_dict()
      new_dict = style.to_dict()

      for key, new_value in new_dict.items():
        if key not in ['query','ref_id']:
          prev_value = prev_dict.get(key)
          if new_value != prev_value:
            handler_func_name = f"transition_{key}"
            handler = getattr(self,handler_func_name)
            if handler:
              handler(ref,style)



class MapStyleController(Controller):
  """
  Manage map styles. Tightly coupled to viewer.
  Expects viewer to be parent
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.state.signals.style_change.connect(self.apply_from_json)

  
  def apply_from_json(self,json_str):
    style = Style.from_json(json_str)
    self.apply(style)

  def apply(self, style):
    # Apply a style

    # Get the ref the style applies to
    ref = self.state.references[style.ref_id]
    if isinstance(ref,MapRef):
      prev_dict = ref.style.to_dict() # old style dict
      new_dict = style.to_dict() # new style dict

      for key, new_value in new_dict.items():
        prev_value = prev_dict.get(key)
        if new_value != prev_value:
          handler_func_name = f"transition_{key}"

          handler = getattr(self,handler_func_name)
          if handler:
            handler(ref,style)


  def transition_iso(self,ref,style,queue=False):
    self.parent.set_iso(ref,style.iso)