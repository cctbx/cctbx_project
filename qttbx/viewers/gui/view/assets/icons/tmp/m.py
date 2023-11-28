from PIL import Image, ImageDraw

def draw_stylized_mandelbrot(draw, size):
  # Define the Mandelbrot outline
  mandelbrot_outline = [
    (20, 20),
    (70, 20),
    (80, 35),
    (60, 50),
    (80, 65),
    (70, 80),
    (20, 80),
    (20, 20),
  ]

  # Define the 'bulb' on the left
  bulb_outline = [
    (30, 40),
    (20, 50),
    (30, 60),
    (40, 50),
    (30, 40),
  ]

  # Draw the Mandelbrot outline
  draw.polygon(mandelbrot_outline, outline="black", fill="black")

  # Draw the 'bulb' on the left
  draw.polygon(bulb_outline, outline="white", fill="white")

# Create an image
size = (100, 100)
image = Image.new("RGB", size, "white")
draw = ImageDraw.Draw(image)

# Draw stylized Mandelbrot
draw_stylized_mandelbrot(draw, size)

# Show image
image.show()

