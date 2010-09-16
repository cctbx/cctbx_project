  bool verbose = false;

  static std::size_t ok_counter = 0;
  static std::size_t error_counter = 0;

  template <typename VectorType1, typename VectorType2>
  void verify(long line, VectorType1 const& v, VectorType2 const& a)
  {
    if (v.size() != a.size()) {
      std::cout << line << ": size mismatch: "
                << v.size() << ", " << a.size() << std::endl;
      error_counter++;
      return;
    }
    else {
      for(std::size_t i=0;i<v.size();i++) {
        if (v[i] != a[i]) {
          std::cout << line << ": value mismatch, index " << i << std::endl;
          error_counter++;
          return;
        }
      }
    }
    if (verbose) std::cout << line << ": OK" << std::endl;
    ok_counter++;
  }

  void check_true(long line, bool stat)
  {
    if (!stat) {
      std::cout << line << ": Error" << std::endl;
      error_counter++;
    }
    else {
      if (verbose) std::cout << line << ": OK" << std::endl;
      ok_counter++;
    }
  }

  void check_false(long line, bool stat) {
    check_true(line, !stat);
  }

  bool approx_equal(double x, double y, double tolerance=1.e-5)
  {
    x -= y;
    if (x < 0) x = -x;
    if (x <= tolerance) return true;
    return false;
  }

#pragma clang diagnostic ignored "-Wunused-function"
