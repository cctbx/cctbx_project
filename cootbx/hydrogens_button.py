from __future__ import absolute_import, division, print_function

class hydrogen_toggle(object):
  def __init__(self, separator=False):
    import coot # import dependency
    import coot_python
    import gtk
    toolbar = coot_python.main_toolbar()
    assert (toolbar is not None)
    if (separator):
      toolbar.insert(gtk.SeparatorToolItem(), -1)
    self.h_button = gtk.ToggleToolButton()
    self.h_button.set_label("Hydrogens off")
    self.h_button.set_is_important(True)
    toolbar.insert(self.h_button, -1)
    self.h_button.connect("clicked", self.OnToggleHydrogens)
    self.h_button.set_active(True)
    self.h_button.show()

  def OnToggleHydrogens(self, *args):
    import coot # import dependency
    if self.h_button.get_active():
      self.h_button.set_label("Hydrogens on")
      for imol in model_molecule_list():
        set_draw_hydrogens(imol, True)
    else :
      self.h_button.set_label("Hydrogens off")
      for imol in model_molecule_list():
        set_draw_hydrogens(imol, False)

if (__name__ == "__main__"):
  hydrogen_toggle(separator=True)
