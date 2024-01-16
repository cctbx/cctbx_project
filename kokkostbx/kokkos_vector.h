#ifndef KOKKOSTBX_VECTOR_H
#define KOKKOSTBX_VECTOR_H

#include <Kokkos_Core.hpp>
#include <iostream>
#include <type_traits>

// Bandaid since parts of diffBragg don't support Kokkos yet
// #ifndef KOKKOS_FUNCTION
// #define KOKKOS_FUNCTION
// #endif

// #ifdef KOKKOS_CORE_HPP
//     template <typename T> T KOKKOS_FUNCTION sqrt_func(T x) { return
//     ::Kokkos::Experimental::sqrt(x); }
// #else
//     #include <cmath>
//     template <typename T> T KOKKOS_FUNCTION sqrt_func(T x) { return sqrt(x); }
// #endif

namespace {
template <typename T>
KOKKOS_FUNCTION typename std::enable_if<std::is_integral<T>::value, void>::type print_num(
    const T& x) {
    printf("%d ", x);
}

template <typename T>
KOKKOS_FUNCTION typename std::enable_if<std::is_floating_point<T>::value, void>::type print_num(
    const T& x) {
    printf("%f ", x);
}
}  // namespace

namespace kokkostbx {

template<typename Derived, typename NumType>
struct stream_initializer {
    Derived& vector;
    int pointer;

    KOKKOS_INLINE_FUNCTION stream_initializer(Derived& v, NumType val) : vector(v), pointer(0) {
        vector[pointer] = val;
    }

    KOKKOS_INLINE_FUNCTION stream_initializer(stream_initializer&& other) = default;

    KOKKOS_FUNCTION stream_initializer& operator,(NumType val) {
        assert(pointer<(int)vector.get_size() &&
               "Too many coefficients to initialize with '<<'-operator'!");
        pointer += 1;
        vector[pointer] = val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION ~stream_initializer() {
        assert((pointer+1)==vector.get_size() &&
                "Too few coefficients to initialize with '<<'-operator!");
    }
};

template <typename Derived, typename NumType, size_t size>
struct vector_base {
    NumType data[size] = {};

    // CONSTRUCTOR
    vector_base() = default;

    KOKKOS_FUNCTION vector_base(NumType val) {
        for (NumType& d : data) {
            d = val;
        }
    }

    KOKKOS_FUNCTION vector_base(NumType arr[]) {
        for (size_t i = 0; i < size; ++i) {
            data[i] = arr[i];
        }
    }

    KOKKOS_FUNCTION vector_base(Derived v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] = v[i];
        }
    }

    // OPERATORS
    // streaming
    friend std::ostream& operator<<(
        std::ostream& os,
        const vector_base<Derived, NumType, size>& v) {
        for (size_t i = 0; i < size; ++i) {
            if (i > 0) {
                os << " ";
            }
            os << v.data[i];
        }
        return os;
    }

    // access
    friend KOKKOS_FUNCTION stream_initializer<Derived, NumType> operator<<(Derived& v, NumType val) {
        return stream_initializer<Derived, NumType>(v, val);
    }


    KOKKOS_FUNCTION NumType& operator[](const int index) { return data[index]; }

    KOKKOS_FUNCTION NumType operator[](const int index) const { return data[index]; }

    // addition
    KOKKOS_FUNCTION friend Derived operator+(const Derived& lhs, const Derived& rhs) {
        Derived sum = lhs;
        sum += rhs;
        return sum;
    }

    // KOKKOS_FUNCTION friend Derived operator+(const Derived& lhs, NumType& rhs) {
    //     Derived sum = lhs;
    //     sum += rhs;
    //     return sum;
    // }

    // KOKKOS_FUNCTION friend Derived operator+(NumType lhs, const Derived& rhs) {
    //     return rhs + lhs;
    // }

    KOKKOS_FUNCTION void operator+=(const Derived& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] += v.data[i];
        }
    }

    KOKKOS_FUNCTION void operator+=(const NumType& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] += v;
        }
    }

    // subtraction
    KOKKOS_FUNCTION friend Derived operator-(const Derived& vec) {
        Derived negative = vec;
        for (size_t i = 0; i < size; ++i) {
            negative[i] *= -1;
        }
        return negative;
    }

    KOKKOS_FUNCTION friend Derived operator-(const Derived& lhs, const Derived& rhs) {
        Derived sum = lhs;
        sum -= rhs;
        return sum;
    }

    KOKKOS_FUNCTION void operator-=(const Derived& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] -= v.data[i];
        }
    }

    KOKKOS_FUNCTION void operator-=(const NumType& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] -= v;
        }
    }

    // multiplication
    KOKKOS_FUNCTION friend Derived operator*(const Derived& lhs, const Derived& rhs) {
        Derived prod = lhs;
        prod *= rhs;
        return prod;
    }

    // KOKKOS_FUNCTION friend Derived operator*(const Derived& lhs, NumType& rhs) {
    //     Derived prod = lhs;
    //     prod *= rhs;
    //     return prod;
    // }

    // KOKKOS_FUNCTION friend Derived operator*(NumType lhs, const Derived& rhs) {
    //     return rhs * lhs;
    // }

    KOKKOS_FUNCTION void operator*=(const Derived& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] *= v.data[i];
        }
    }

    KOKKOS_FUNCTION void operator*=(const NumType& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] *= v;
        }
    }

    // division
    KOKKOS_FUNCTION friend Derived operator/(const Derived& lhs, const Derived& rhs) {
        Derived quot = lhs;
        quot /= rhs;
        return quot;
    }

    KOKKOS_FUNCTION void operator/=(const Derived& v) {
        for (size_t i = 0; i < size; ++i) {
            data[i] /= v.data[i];
        }
    }

    KOKKOS_FUNCTION void operator/=(const NumType& v) {
        const NumType v_r = 1 / v;
        for (size_t i = 0; i < size; ++i) {
            data[i] *= v_r;
        }
    }

    // METHODS
    KOKKOS_FUNCTION void print(const char name[]) const {
        printf("%s: ", name);
        for (size_t i = 0; i < size; ++i) {
            print_num(data[i]);
        }
        printf("\n");
    }

    static constexpr KOKKOS_FUNCTION size_t get_size() { return size; }

    KOKKOS_FUNCTION void zero() {
        for (size_t i = 0; i < size; ++i) {
            data[i] = 0;
        }
    }

    KOKKOS_FUNCTION void ones() {
        for (size_t i = 0; i < size; ++i) {
            data[i] = 1;
        }
    }

    KOKKOS_FUNCTION bool is_zero() const {
        for (size_t i = 0; i < size; ++i) {
            if (data[i] != 0)
                return false;
        }
        return true;
    }

    KOKKOS_FUNCTION NumType dot(const Derived& v) const {
        NumType sum = 0;
        for (size_t i = 0; i < size; ++i) {
            sum += data[i] * v.data[i];
        }
        return sum;
    }

    KOKKOS_FUNCTION NumType length_sqr() const {
        NumType sum = 0;
        for (size_t i = 0; i < size; ++i) {
            sum += data[i] * data[i];
        }
        return sum;
    }

    KOKKOS_FUNCTION NumType length() const {
        // return sqrt_func(length_sqr());
        return ::Kokkos::sqrt(length_sqr());
    }

    KOKKOS_FUNCTION void normalize() {
        NumType l = length();
        if (l > 0) {
            NumType l_r = 1 / l;
            for (size_t i = 0; i < size; ++i) {
                data[i] *= l_r;
            }
        }
    }

    KOKKOS_FUNCTION Derived get_unit_vector() const {
        NumType l = length();
        Derived unit_vector{};
        if (l > 0) {
            NumType l_r = 1 / l;
            for (size_t i = 0; i < size; ++i) {
                unit_vector[i] = data[i] * l_r;
            }
        }
        return unit_vector;
    }
};

template <typename NumType, size_t size>
struct vector : public vector_base<vector<NumType, size>, NumType, size> {
    using vector_base = kokkostbx::vector_base<vector<NumType, size>, NumType, size>;

    vector() = default;
    KOKKOS_FUNCTION vector(NumType val) : vector_base(val){};
    KOKKOS_FUNCTION vector(NumType arr[]) : vector_base(arr){};
};

}  // namespace kokkostbx

#endif
