#ifndef IOTBX_CIF_LOOP_FORMAT_H
#define IOTBX_CIF_LOOP_FORMAT_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <string>

namespace iotbx { namespace cif {

  namespace af = scitbx::af;

  //! Rows start to stop of the columns as iotbx.cif.model.loop.show()
  //! prints them, each followed by a newline.
  /*! columns hold the printed values (quoting applied) and have the same
      size.  widths has one entry per column, the field width as the %Ns
      of a printf format: positive right-aligns, negative left-aligns,
      zero pads nothing; padding counts UTF-8 code points.  The padded
      values follow the prefix joined by two spaces, and the line loses
      its trailing whitespace as str.rstrip() removes it.  An empty widths
      selects the non-aligned form: prefix, the values joined by one
      space, no strip.  std::invalid_argument on a bad range, columns of
      unequal size, or a widths size other than the column count.
   */
  std::string
  format_loop_rows(
    af::const_ref<af::shared<std::string> > const& columns,
    std::size_t start,
    std::size_t stop,
    std::string const& prefix,
    af::const_ref<int> const& widths);

  //! True if a value of the column holds a newline or is changed by
  //! iotbx.cif.model.format_value, whose rules this ports: change both
  //! together (exercise_quoting_scan).
  bool
  column_needs_format_value(af::const_ref<std::string> const& column);

  //! Index of the first value at or after start that is neither ".", "?"
  //! nor a plain number; the size of the column when there is none.
  /*! A plain number is an optional sign, digits with an optional fraction
      or a fraction alone, and an optional exponent of at most two digits,
      at most 40 characters in all: flex.double(flex.std_string) accepts
      every such value and none of them overflows or underflows.
   */
  std::size_t
  first_non_plain_number(
    af::const_ref<std::string> const& column,
    std::size_t start);

}} // namespace iotbx::cif

#endif // IOTBX_CIF_LOOP_FORMAT_H
