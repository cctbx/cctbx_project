
import matplotlib.colors as mcolors

class Color:
  def __init__(self, color):
    if isinstance(color, str):
      self.from_string(color)
    elif isinstance(color, tuple) and len(color) == 3:
      if all(isinstance(c, int) for c in color):
        self.from_rgb_int(color)
      elif all(isinstance(c, float) for c in color):
        self.from_rgb_float(color)
    else:
      raise ValueError("Invalid color format")

  def from_string(self, color):
    if color in mcolors.CSS4_COLORS:
      self.rgb = mcolors.to_rgb(mcolors.CSS4_COLORS[color])
    elif color.startswith('#'):
      self.rgb = mcolors.to_rgb(color)
    else:
      raise ValueError("Unknown named color or invalid hex code")

  def from_rgb_int(self, rgb):
    if all(0 <= c <= 255 for c in rgb):
      self.rgb = tuple(c / 255 for c in rgb)
    else:
      raise ValueError("RGB integer values must be between 0 and 255")

  def from_rgb_float(self, rgb):
    if all(0.0 <= c <= 1.0 for c in rgb):
      self.rgb = rgb
    else:
      raise ValueError("RGB float values must be between 0.0 and 1.0")

  def __repr__(self):
    return f"Color(rgb={self.rgb})"

  @property
  def r(self):
    return self.rgb[0]

  @property
  def g(self):
    return self.rgb[1]

  @property
  def b(self):
    return self.rgb[2]

  @property
  def R(self):
    return int(self.rgb[0] * 255)

  @property
  def G(self):
    return int(self.rgb[1] * 255)

  @property
  def B(self):
    return int(self.rgb[2] * 255)

  @property
  def RGB(self):
    return (self.R, self.G, self.B)