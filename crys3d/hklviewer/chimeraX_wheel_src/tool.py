from __future__ import division

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.tools import ToolInstance
import math

class HKLviewerTool(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.
    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:user/tools/tutorial.html"
                                # Let ChimeraX know about our help page


    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "HKLviewer"
        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()
        self.clipper_crystdict = None
        session.HKLviewer = self

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QHBoxLayout
        from . import HKLviewer

        hbox = QHBoxLayout()
        self.Guiobj = HKLviewer.run(isembedded=True, chimeraxsession=self.session)
        hbox.addWidget(self.Guiobj.window)
        self.tool_window.ui_area.setLayout(hbox)
        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def isolde_clipper_data_to_dict(self):
        try:
            sh = self.session.isolde.selected_model.parent
            xmapset = sh.map_mgr.xmapsets[0]
            clipperlabel= list(xmapset.experimental_data.items())[0][0]
            labels = [lbl.strip() for lbl in clipperlabel.split(",")]
            self.clipper_crystdict = {}
            self.clipper_crystdict["spg_number"] = xmapset.spacegroup.spacegroup_number
            self.clipper_crystdict["unit_cell"] =  (xmapset.unit_cell.cell.a, xmapset.unit_cell.cell.b, xmapset.unit_cell.cell.c,
               xmapset.unit_cell.cell.alpha*180/math.pi, xmapset.unit_cell.cell.beta*180/math.pi, xmapset.unit_cell.cell.gamma*180/math.pi)
            self.clipper_crystdict["HKL"] = xmapset.experimental_data[clipperlabel].data.data[0].tolist()
            self.clipper_crystdict[clipperlabel] = xmapset.experimental_data[clipperlabel].data.data[1].transpose().tolist()
            self.clipper_crystdict["FCALC,PHFCALC"] = xmapset.live_xmap_mgr.f_calc.data[1].transpose().tolist()
            self.clipper_crystdict["2FOFC,PH2FOFC"] = xmapset.live_xmap_mgr.base_2fofc.data[1].transpose().tolist()
        except Exception as e:
            pass


    def return_pressed(self):
        # The use has pressed the Return key; log the current text as HTML
        from chimerax.core.commands import run
        # ToolInstance has a 'session' attribute...
        run(self.session, "log html %s" % self.line_edit.text())

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.)
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from PyQt5.QtWidgets import QAction
        settings_action = QAction("HKLviewer settings", menu)
        settings_action.triggered.connect(lambda *args: self.Guiobj.SettingsDialog() )
        menu.addAction(settings_action)

    def delete(self):
      from Qt.QtCore import QEvent
      self.session.triggers.remove_handler(self.Guiobj.chimeraxprocmsghandler)
      self.Guiobj.closeEvent(QEvent.Close)
      super(HKLviewerTool, self).delete()

    def take_snapshot(self, session, flags):
        return {
            'version': 1,
            'current text': self.line_edit.text()
        }

    @classmethod
    def restore_snapshot(class_obj, session, data):
        # Instead of using a fixed string when calling the constructor below, we could
        # have saved the tool name during take_snapshot() (from self.tool_name, inherited
        # from ToolInstance) and used that saved tool name.  There are pros and cons to
        # both approaches.
        inst = class_obj(session, "cctbx.HKLviewer")
        inst.line_edit.setText(data['current text'])
        return inst
