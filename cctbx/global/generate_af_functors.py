def generate_unary_functor(name, op):
  print """
  template <typename ResultType,
            typename ArgumentType>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType const& x) const {
      return ResultType(%s);
    }
  };""" % (name, op)

def generate_binary_functor(name, op):
  print """
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y) const {
      return ResultType(%s);
    }
  };""" % (name, op)

def generate_in_place_binary_functor(name, op):
  print """
  template <typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    ArgumentType1& operator()(ArgumentType1& x,
                              ArgumentType2 const& y) const {
      %s;
      return x;
    }
  };""" % (name, op)

def generate_greater_less_functor(name, op):
  print """
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y) const {
      return ResultType(x %s y);
    }
    ResultType reverse(ArgumentType2 const& y,
                       ArgumentType1 const& x) const {
      return ResultType(y %s x);
    }
  };""" % (name, op, op)
