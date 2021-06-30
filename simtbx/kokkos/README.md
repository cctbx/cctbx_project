# simtbx - Kokkos
*Note to developers*

It is problematic to capture a class member using a lambda function in Kokkos. For more details, see <https://github.com/kokkos/kokkos/issues/695>.

The basic problem is that `m_member` gets turned by the lambda into `this->m_member` and then the host pointer is copied to the device, where it raises a `cudaErrorIllegalAddress` at runtime.

There are a number of options around this, we use option 1 or 3:

1. Create a local variable `auto local_member = m_member;` and only use these in a lambda 
2. Use C++17 and the `KOKKOS_CLASS_LAMBDA` instead of `KOKKOS_LAMBDA`. However, this copies the whole class to the device and fails if any member has no device equivalent, like `std::vector`
3. Don't use member functions, use outside functions
4. Don't use lambdas, use functors

More information: <https://github.com/kokkos/kokkos/wiki/Lambda-Dispatch#lambdas-inside-classes>
