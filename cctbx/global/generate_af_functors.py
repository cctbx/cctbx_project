def generate_unary_functor(name, op):
  print """
  template <typename ResultType,
            typename ArgumentType>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(const ArgumentType& x) const {
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
    ResultType operator()(const ArgumentType1& x,
                          const ArgumentType2& y) const {
      return ResultType(%s);
    }
  };""" % (name, op)

def generate_in_place_binary_functor(name, op):
  print """
  template <typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    ArgumentType1& operator()(ArgumentType1& x,
                              const ArgumentType2& y) const {
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
    ResultType operator()(const ArgumentType1& x,
                          const ArgumentType2& y) const {
      return ResultType(x %s y);
    }
    ResultType reverse(const ArgumentType2& y,
                       const ArgumentType1& x) const {
      return ResultType(y %s x);
    }
  };""" % (name, op, op)
