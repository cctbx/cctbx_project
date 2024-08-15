"""
A common place to define and convert between colors among many graphics programs.
"""
import matplotlib.colors as mcolors

class Color:
  """
  A single color. Various methods to construct:
    1. rgb float
    2. rgb int
    3. html string
    4. label string
    5. transparency float
  """
  def __init__(self, r,g,b,transparency=1.0):
    # a tuple of floats (r,g,b) is the internal state
    rgb = (r,g,b)
    if all(0.0 <= c <= 1.0 for c in rgb):
      self.rgb = rgb
      self.transparency = transparency
    else:
      raise ValueError("RGB float values must be between 0.0 and 1.0")

  @classmethod
  def from_string(cls, color,transparency=1.0):
    if color in mcolors.CSS4_COLORS:
      rgb = mcolors.to_rgb(mcolors.CSS4_COLORS[color])
      return cls(*rgb,transparency=transparency)
    elif color.startswith('#'):
      rgb = mcolors.to_rgb(color)
      return cls(*rgb,transparency=transparency)
    else:
      raise ValueError("Unknown named color or invalid hex code")

  @classmethod
  def from_rgb_int(cls, rgb,transparency=255):
    if all(0 <= c <= 255 for c in (*rgb,transparency)):
      # Convert to float for single instantiation
      rgb = tuple(c / 255 for c in rgb)
      transparency = transparency/255
      return cls(*rgb,transparency=transparency)
    else:
      raise ValueError("RGB integer values must be between 0 and 255")

  @classmethod
  def from_rgb_float(cls,rgb,transparency=1.0):
    return cls(rgb,transparency=transparency)
  
  @classmethod
  def from_rgb(cls,rgb,transparency=1.0):
    if isinstance(rgb[0],int):
      return cls.from_rgb_int(rgb,int(transparency*255))
    else:
      return cls.from_rgb_float(rgb,transparency=transparency)

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
  def a(self):
    return self.transparency

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
  def A(self):
    return int(self.transparency * 255)

  @property
  def RGB(self):
    return (self.R, self.G, self.B)

  @property
  def RGBA(self):
    return (self.R, self.G, self.B, self.A)