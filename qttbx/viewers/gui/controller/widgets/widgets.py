from PySide2.QtWidgets import QMessageBox


from ..controller import Controller
from dataclasses import replace



class ISOWidgetController(Controller):
  def __init__(self,parent=None,view=None,map_ref=None):
    super().__init__(parent=parent,view=view)
    self.map_ref = map_ref
    self.view.map_ref = map_ref

    # Connect
    self.view.slider.valueChanged.connect(self.iso_update)


  def iso_update(self,value):
    #print("slider update: active_map_ref",self.map_ref.id)
    real_value = self.view.slider.slider_iso_to_real(value,self.view.slider.slider_iso_granularity)
    self.view.parent().iso_label.setText(f"ISO: {round(real_value,2)}")
    style = replace(self.map_ref.style,iso=real_value)
    self.map_ref.style = style




class OpacityWidgetController(Controller):
  def __init__(self,parent=None,view=None,map_ref=None):
    super().__init__(parent=parent,view=view)
    self.map_ref = map_ref
    self.view.slider.valueChanged.connect(self.opacity_update)


  def opacity_update(self,value):
    ref_id = self.map_ref.id
    print("slider update opacity:",ref_id,value)
    QMessageBox.information(self.view, 'Notice', "Opacity slider not yet connected")
    real_value = value/100
    self.view.label.setText(f"Opacity: {round(real_value,2)}")
