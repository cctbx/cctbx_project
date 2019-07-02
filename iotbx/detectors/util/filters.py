from __future__ import absolute_import, division, print_function
from six.moves import range
import math,copy
from scitbx.array_family import flex
from iotbx.detectors import low_pass_filter
from libtbx import adopt_init_args

def get_factorizable_block(asic):
  input1 = asic[2]-asic[0]
  input2 = asic[3]-asic[1]
  return (asic[0],asic[1],asic[0]+primecheck_minus(input1),
                          asic[1]+primecheck_minus(input2))

def primecheck_minus(input_m):
  residual = copy.copy(input_m)
  return_m = input_m
  while (residual%2==0): residual//=2;
  while (residual%3==0): residual//=3;
  while (residual%5==0): residual//=5;
  if (residual>1):
    return_m = primecheck_minus(return_m-1);
  return return_m;

def primecheck_plus(input_m):
  residual = copy.copy(input_m)
  return_m = input_m
  while (residual%2==0): residual//=2;
  while (residual%3==0): residual//=3;
  while (residual%5==0): residual//=5;
  if (residual>1):
    return_m = primecheck_plus(return_m+1);
  return return_m;

def hi_pass_filter(complex_data):

  sz_x = complex_data.focus()[0]
  sz_y = complex_data.focus()[1]
  for x in range(sz_x):
    #dx = min(x,sz_x-x)
    dx = abs((sz_x/2)-x)
    for y in range(sz_y):
      #dy = min(y,sz_y-y)
      dy = abs((sz_y/2)-y)
      dist = math.hypot(dx,dy)
      if dist<100.:  complex_data[(x,y)]=complex(0.,0.)
      #else:
      #  complex_data[(x,y)]*=math.exp(-dist*dist/20.)

class padded_unpadded:
  def __init__(self, data, asic):
    adopt_init_args(self,locals())

    padding = 10 # minimum padding pixels on each side
    self.input1 = input1 = asic[2]-asic[0]
    self.input2 = input2 = asic[3]-asic[1]

    self.size1 = size1 = primecheck_plus(input1+2*padding)
    self.leading1 = (size1-input1)//2
    self.trailing1 = size1-self.leading1-input1

    self.size2 = size2 = primecheck_plus(input2+2*padding)
    self.leading2 = (size2-input2)//2
    self.trailing2 = size2-self.leading2-input2

  def get_padded_input_data(self):
    result = self.data.__class__(flex.grid(self.size1,self.size2))

    result.matrix_paste_block_in_place(
      block = self.data.matrix_copy_block(
        i_row=self.asic[0],i_column=self.asic[1],
        n_rows=self.input1,
        n_columns=self.input2),
      i_row = self.leading1,
      i_column = self.leading2
    )

    leading_line = result.matrix_copy_block(
      i_row=self.leading1,i_column=0,n_rows=1,n_columns=self.size2)
    for irow in range(self.leading1):
      result.matrix_paste_block_in_place(
        block = leading_line,i_row=irow,i_column=0)

    trailing_line = result.matrix_copy_block(
      i_row=self.size1-self.trailing1-1,i_column=0,n_rows=1,n_columns=self.size2)
    for irow in range(self.trailing1):
      result.matrix_paste_block_in_place(
        block = trailing_line,i_row=self.size1-self.trailing1+irow,i_column=0)

    leading_line = result.matrix_copy_block(
      i_row=0,i_column=self.leading2,n_rows=self.size1,n_columns=1)
    for icolumn in range(self.leading2):
      result.matrix_paste_block_in_place(
        block = leading_line,i_row=0,i_column=icolumn)

    trailing_line = result.matrix_copy_block(
      i_row=0, i_column=self.size2-self.trailing2-1,n_rows=self.size2,n_columns=1)
    for icolumn in range(self.trailing2):
      result.matrix_paste_block_in_place(
        block = trailing_line,i_row=0, i_column=self.size2-self.trailing2+icolumn)

    return result

  def get_unpadded_result_data(self,corrected_data):
    return corrected_data.matrix_copy_block(
        i_row=self.leading1,i_column=self.leading2,
        n_rows=self.input1,
        n_columns=self.input2)

def background_correct_padded_block(data, raw_asic):

  Pad = padded_unpadded(data,raw_asic)

  block = Pad.get_padded_input_data()

  complex_data = flex.polar(block.as_double(),flex.double(flex.grid(block.focus())))
  from scitbx import fftpack
  fft = fftpack.complex_to_complex_2d(block.focus())
  # input data here
  fft.forward(complex_data)
  # your manipulation here
  low_pass_filter(complex_data)

  fft.backward(complex_data)
  # real data
  filtered_data = flex.real(complex_data)/(fft.n()[0]*fft.n()[1])
  # XXX change this depending on float/int data type:
  corrected_data = block - filtered_data.iround()

  return Pad.get_unpadded_result_data(corrected_data)

def background_correct(data, raw_asic):

  prime_asic = get_factorizable_block(raw_asic)
  print("Working on block",prime_asic)
  block = data.matrix_copy_block(
      i_row=prime_asic[0],i_column=prime_asic[1],
      n_rows=prime_asic[2]-prime_asic[0],
      n_columns=prime_asic[3]-prime_asic[1])

  complex_data = flex.polar(block.as_double(),flex.double(flex.grid(block.focus())))
  from scitbx import fftpack
  fft = fftpack.complex_to_complex_2d(block.focus())
  # input data here
  fft.forward(complex_data)
  # your manipulation here
  low_pass_filter(complex_data)

  fft.backward(complex_data)
  # real data
  filtered_data = flex.real(complex_data)/(fft.n()[0]*fft.n()[1])
  corrected_data = block - filtered_data.iround()

  requested_block = data.matrix_copy_block(
      i_row=raw_asic[0],i_column=raw_asic[1],
      n_rows=raw_asic[2]-raw_asic[0],
      n_columns=raw_asic[3]-raw_asic[1])
  requested_block.matrix_paste_block_in_place(
      block = corrected_data,
      i_row = prime_asic[0] - raw_asic[0],
      i_column = prime_asic[1] - raw_asic[1]
      )

  return requested_block
