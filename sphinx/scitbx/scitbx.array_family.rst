========================
scitbx.array_family.flex
========================

The array family lies at the core of all numerical routines in CCTBX, and is
designed to make the transition between Python and C++ code as seamless as
possible.  The Python API mimics that of the built-in list type, with many
extensions specific to the types involved.  At the core of the array family is
the :py:mod:`~scitbx.array_family.flex` submodule, which includes most of the array types, including:

* :py:class:`.flex.int` - maps to C++ type :cpp:type:`int`
* :py:class:`.flex.long` - maps to C++ type :cpp:type:`long`
* :py:class:`.flex.size_t` - maps to C++ type :cpp:type:`std::size_t`
* :py:class:`.flex.bool` - maps to C++ type :cpp:type:`bool`
* :py:class:`.flex.double` - maps to C++ type :cpp:type:`double`
* :py:class:`.flex.float` - maps to C++ type :cpp:type:`float`
* :py:class:`.flex.complex_double` - maps to C++ type :cpp:type:`std::complex<double>`
* :py:class:`.flex.std_string` - maps to C++ type :cpp:type:`std::string`
* :py:class:`.flex.vec2_double` - An array of 2-tuple doubles (e.g. for 2D coordinates)
* :py:class:`.flex.vec3_double` - An array of 3-tuple doubles (e.g. for 3D coordinates)
* :py:class:`.flex.vec3_int` - An array of 3-tuple ints (e.g. for 3D grid coordinates)
* :py:class:`.flex.sym_mat3_double` - An array of 6-tuple doubles. Used to store arrays
  of symmetric 3x3 matrices

>>> from scitbx.array_family import flex
>>> I = flex.int([1,2,3,4])
>>> I.append(5)
>>> I.extend(flex.int([6]))
>>> len(I)
6
>>> del I[0]
>>> print list(I)
[2, 3, 4, 5, 6]
>>> print list(I[2:4])
[4, 5]
>>> I.all_ne(1)
True
>>> print list(I.as_double())
[2.0, 3.0, 4.0, 5.0, 6.0]
>>> J = flex.int(5, 4)
>>> list(J)
[4, 4, 4, 4, 4]

For numeric element types the flex type supports algebraic operations::

  >>> a = flex.double([1,2,3])
  >>> b = flex.double([3,2,1])
  >>> tuple(a+b)
  (4.0, 4.0, 4.0)
  >>> tuple(flex.sqrt(a+b))
  (2.0, 2.0, 2.0)

The flex type also supports a large number of other functions (abs, sin, pow, etc.),
slicing, and importantly, pickling_.

.. _pickling https://docs.python.org/2.7/library/pickle.html

If all this looks similar to the popular NumPy package: it is at the surface.
However, there are two very important differences:

1. Under the hood the flex types are instantiations of a C++ array type that
resembles familiar STL container types as much as possible. In contrast
NumPy presents itself with a raw 'C' API.

2. It is straightforward to implement other flex types with custom user-defined
element types, even outside the scitbx module. This would be extremely
difficult to do with NumPy, and is virtually impossible if the user-defined
types are implemented in C++.

The extensible nature of the :py:module:`scitbx.array_family` library is used
elsewhere in CCTBX to wrap everything from additional
crystallography-specific numerical types such as Miller indices (``(h,k,l)``),
to X-ray scatterer objects, atoms represented in a PDB file, and geometry
restraints.  These array types are described elsewhere, but many of the same
principles explained here (especially selections) will apply to them as well.

Create a toy array::

  a = flex.double([10,11,12])

One way of looping over the array::

  for i in range(a.size()):
    print a[i]

A better way of looping over the array::

  for ai in a:
    print ai

Another good way of looping over the array::

  for i,ai in enumerate(a):
    print i, ai

Modify the elements one-by-one::

  for i in range(a.size()):
    a[i] *= 10
  for ai in a:
    print ai

A better way of modifying all elements::

  a += 100 # this works at C++ speed
  for ai in a:
    print ai


Multidimensional arrays: flex.grid
----------------------------------

The underlying data stored in a flex array is always stored in a contiguous
one-dimensional sequence of memory. Multidimensional flex arrays
(up to 10 dimensions are supported) are provided by attaching an accessor,
of type :py:class:`.flex.grid` to a flex array, which specifies how the
underlying one-dimensional array of data should be interpreted as a
multidimensional array.

::

  >>> a = flex.double(range(9))
  >>> grid = a.accessor()
  >>> print grid.nd() # the number of dimensions
  1
  >>> print grid.all() # the size in each dimension
  (9,)

Resize an existing array::

  >>> a.reshape(flex.grid(3,3))
  >>> print a.nd()
  2
  >>> print a.all()
  (3,3)

Copy an accessor from one array to another::

  >>> b = flex.int(range(9))
  >>> print b.nd()
  1
  >>> b.reshape(a.accessor())
  >>> print b.nd()
  2

Multidimensional flex arrays are stored in `row-major order`_ and elements
can be accessed using the Python square bracket notation::

  >>> c = flex.int(range(6))
  >>> c.reshape(flex.grid(2,3))
  >>> print list(c)
  [0, 1, 2, 3, 4, 5]
  >>> for i in range(c.all()[0]):
  ...   for j in range(c.all()[1]):
  ...     print c[i,j],
  ...   print
  ...
  0 1 2
  3 4 5

.. _row-major order: http://en.wikipedia.org/wiki/Row-major_order

Python's slice_ notation also works for multidimensional arrays::

  >>> print list(c[0:2, 0:2])
  [0, 1, 3, 4]
  >>> c[0:2, 0:2] += 10
  >>> print list(c)
  [10, 11, 2, 13, 14, 5]

.. _slice: https://docs.python.org/2/library/functions.html?highlight=slice#slice


Working with selections
-----------------------

A particularly powerful feature of flex arrays is the ability to select a subset
of array values, using either :py:class:`.flex.bool` or :py:class:`.flex.size_t`
to indicate the elements desired::

  >>> from scitbx.array_family import flex
  >>> I = flex.int(range(10,30))
  >>> sel = flex.bool()
  >>> for x in range(20) :
  ...   if (x \% 4 == 0) :
  ...     sel.append(True)
  ...   else :
  ...     sel.append(False)
  ...
  >>> J = I.select(sel)
  >>> type(J)
  <class 'scitbx_array_family_flex_ext.int'>
  >>> list(J)
  [10, 14, 18, 22, 26]
  >>> inverse_sel = ~sel
  >>> inverse_isel = inverse_sel.iselection()
  >>> type(inverse_isel)
  <class 'scitbx_array_family_flex_ext.size_t'>
  >>> list(inverse_isel)
  [1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19]
  >>> K = I.select(inverse_isel)
  [11, 12, 13, 15, 16, 17, 19, 20, 21, 23, 24, 25, 27, 28, 29]


Here we have generated a boolean selection of the same size as the target
array, used this to pull out a subset of array values (still as a
:py:class:`.flex.int` object), obtained the inverse boolean selection,
converted this to selected array indices instead of boolean flags, and
performed another selection using the indices.  In most situations the
boolean flags and the enumerated indices can be used interchangeably depending
on which is most convenient. For instance, setting selected elements to
desired values::

  >>> I.set_selected(sel, 999)
  <scitbx_array_family_flex_ext.int object at 0x1029eafc8>
  >>> list(I)
  [999, 11, 12, 13, 999, 15, 16, 17, 999, 19, 20, 21, 999, 23, 24, 25, 999, 27, 28, 29]
  >>> print (I.set_selected(inv_isel, -1))
  [999, -1, -1, -1, 999, -1, -1, -1, 999, -1, -1, -1, 999, -1, -1, -1, 999, -1, -1, -1]
  >>> J = flex.int(range(10,20))
  >>> K = flex.int(range(35, 40))
  >>> isel = flex.size_t([5,6,7,8,9])
  >>> J.set_selected(isel, K)
  [10, 11, 12, 13, 14, 35, 36, 37, 38, 39]

(Note that in this example, :code:`array.set_selected()` returns an array -
however, this is the same array modified in place.)

Obviously the selections (and other similar operations) will only work if the
values are compatible; failing to ensure appropriate selections results in
errors::

  >>> from scitbx.array_family import flex
  >>> I = flex.int(range(10))
  >>> sel = flex.bool([ True for x in range(9) ])
  >>> I.select(sel)
  Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
  RuntimeError: scitbx Internal Error:
        /Users/nat/phenix/src/cctbx_project/scitbx/array_family/selections.h(44):
        SCITBX_ASSERT(flags.size() == self.size()) failure.
  >>> isel = flex.size_t([1,5,13])
  >>> I.select(isel)
  Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
  RuntimeError: scitbx Internal Error:
        /Users/nat/phenix/src/cctbx_project/scitbx/array_family/selections.h(21):
        SCITBX_ASSERT(indices[i] < self.size()) failure.

In the first case we have attempted to use an improperly sized :py:class:`flex.bool` array
to denote the selection; in the second, a :py:class:`flex.size_t` array containing
elements with a value that overflows the bounds of the target array.

For selections of type :py:class:`flex.size_t`, the items do not necessarily
have to be in sequential order, and in fact this allows us to reorder any flex
array via a selection::

  >>> from scitbx.array_family import flex
  >>> I = flex.int([10,20,30,40,50,60,70,80,90,100])
  >>> sel = flex.size_t([9,7,5,3,1])
  >>> J = I.select(sel)
  >>> print list(J)
  [100, 80, 60, 40, 20]

We also take advantage of this behavior to enable sorting of flex arrays using
the accessory function :py:class:`flex.sort_permutation`, which generates a
selection rather than directly sorting the array.  Continuing the above
example::

  >>> sel_perm = flex.sort_permutation(J)
  >>> print list(sel_perm)
  [4, 3, 2, 1, 0]
  >>> K = J.select(sel_perm)
  >>> print list(K)
  [20, 40, 60, 80, 100]


API documentation
-----------------

.. automodule:: scitbx.array_family.flex
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.grid
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.int
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.long
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.size_t
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.bool
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.double
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.float
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.complex_double
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.std_string
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.vec2_double
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.vec3_double
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.vec3_int
    :members:
    :undoc-members:

.. autoclass:: scitbx.array_family.flex.sym_mat3_double
    :members:
    :undoc-members:
