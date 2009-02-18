#ifndef BOOST_ADAPTBX_FORWARD_COMPATIBILITY_H
#define BOOST_ADAPTBX_FORWARD_COMPATIBILITY_H

/// Compatibility with older versions of Python


/// Move from int to Py_ssize_t (c.f. PEP 353)
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

#endif // Guard

