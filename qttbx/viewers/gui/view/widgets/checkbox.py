from PySide2.QtWidgets import QCheckBox
from PySide2.QtGui import QMouseEvent

class ConditionalCheckBox(QCheckBox):
    """
    Intercept a check click and only change state upon condition
    """
    def __init__(self, label, parent=None):
        super().__init__(label, parent)

    def mousePressEvent(self, event: QMouseEvent):
        # Condition to check before toggling
        if self.can_toggle():
            super().mousePressEvent(event)
        else:
            # Optionally handle the event if the condition is not met
            #self.log("Condition not met, checkbox state not toggled.")
            pass

    def can_toggle(self):
        # Insert a condition here
        return not self.isChecked()