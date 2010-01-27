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

  template <typename ElementType, typename UnsignedType>
  boost::python::object
  add_selected_unsigned_s(
    boost::python::object const& self,
    af::const_ref<UnsignedType> const& indices,
    ElementType const& value)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < a.size());
      a[indices[i]] += value;
    }
    return self;
  }

  template <typename ElementType>
  boost::python::object
  add_selected_bool_a(
    boost::python::object const& self,
    af::const_ref<bool> const& flags,
    af::const_ref<ElementType> const& values)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    SCITBX_ASSERT(a.size() == flags.size());
    if (values.size() == flags.size()) {
      ElementType* ai = a.begin();
      const bool* fi = flags.begin();
      const ElementType* ni = values.begin();
      const ElementType* ne = values.end();
      while (ni != ne) {
        if (*fi++) *ai += *ni;
        ai++;
        ni++;
      }
    }
    else {
      std::size_t i_value = 0;
      for(std::size_t i=0;i<flags.size();i++) {
        if (flags[i]) {
          SCITBX_ASSERT(i_value < values.size());
          a[i] += values[i_value];
          i_value++;
        }
      }
      SCITBX_ASSERT(i_value == values.size());
    }
    return self;
  }


  template <typename ElementType>
  boost::python::object
  add_selected_bool_s(
    boost::python::object const& self,
    af::const_ref<bool, flex_grid<> > const& flags,
    ElementType const& value)
  {
    boost::python::extract<af::ref<ElementType, flex_grid<> > > a_proxy(self);
    af::ref<ElementType, flex_grid<> > a = a_proxy();
    SCITBX_ASSERT(a.accessor() == flags.accessor());
    for(std::size_t i=0;i<flags.size();i++) {
      if (flags[i]) a[i] += value;
    }
    return self;
  }

}}} // namespace scitbx::af::boost_python
