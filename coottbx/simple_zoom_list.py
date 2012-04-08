
# basic window for displaying a list of features of interest - these could
# be ligands, difference map peaks, validation outliers, etc. (adapted from
# templates in mmtbx.validation scripts).  This will usually be concatenated
# with the actual data of interest, plus whatever additional commands we want
# to run in Coot.

def draw_simple_zoom_list (title, items, zoom_level=30) :
  import gtk
  import coot # import dependency
  window = gtk.Window(gtk.WINDOW_TOPLEVEL)
  scrolled_win = gtk.ScrolledWindow()
  outside_vbox = gtk.VBox(False, 2)
  inside_vbox = gtk.VBox(False, 0)
  window.set_default_size(300, 200)
  window.set_title(title)
  inside_vbox.set_border_width(2)
  window.add(outside_vbox)
  outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
  scrolled_win.add_with_viewport(inside_vbox)
  scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
  frame = gtk.Frame(title)
  vbox = gtk.VBox(False, 0)
  inside_vbox.pack_start(frame, False, False, 2)
  frame.add(vbox)
  def callback_recenter (widget, xyz) :
    set_rotation_centre(*xyz)
    set_zoom(zoom_level)
  for feature, xyz in items :
    button = gtk.Button(feature)
    button.connect("clicked", callback_recenter, xyz)
    vbox.pack_start(button, False, False, 1)
  outside_vbox.set_border_width(2)
  ok_button = gtk.Button("  Close  ")
  outside_vbox.pack_end(ok_button, False, False, 0)
  ok_button.connect("clicked", lambda x: window.destroy())
  window.show_all()
