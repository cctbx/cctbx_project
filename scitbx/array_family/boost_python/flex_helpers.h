namespace scitbx { namespace af { namespace boost_python {

  template <typename ElementType, typename UnsignedType>
  boost::python::object
  add_selected_unsigned_a(
    boost::python::object const& self,
    af::const_ref<UnsignedType> const& indices,
    af::const_ref<ElementType> const& values)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    SCITBX_ASSERT(indices.size() == values.size());
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < a.size());
      a[indices[i]] += values[i];
    }
    return self;
  }

}}} // namespace scitbx::af::boost_python
