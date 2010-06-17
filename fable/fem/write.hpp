#ifndef FEM_WRITE_HPP
#define FEM_WRITE_HPP

#include <fem/common.hpp>
#include <fem/format.hpp>
#include <fem/star.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/real_as_string.hpp>

namespace fem {

  class write_loop // TODO copy-constructor potential performance problem
  {
    private:
      std::auto_ptr<utils::simple_ostream> out;
      int internal_file_len;
      unsigned pos;
      bool prev_was_string;
      io_modes io_mode;
      format::token_loop fmt_loop;
      bool suppress_new_line_at_end;

    public:

      write_loop(
        common& cmn,
        int const& unit,
        unformatted_type const&)
      :
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        pos(0),
        prev_was_string(false),
        io_mode(io_unformatted),
        suppress_new_line_at_end(false)
      {}

      write_loop(
        common& cmn,
        int const& unit,
        star_type const&)
      :
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        pos(0),
        prev_was_string(false),
        io_mode(io_list_directed),
        suppress_new_line_at_end(false)
      {}

      write_loop(
        common& cmn,
        int const& unit,
        char const* fmt)
      :
        out(cmn.io.simple_ostream(unit)),
        internal_file_len(-1),
        pos(0),
        prev_was_string(false),
        io_mode(io_formatted),
        fmt_loop(fmt),
        suppress_new_line_at_end(false)
      {}

      write_loop(
        str_ref const& internal_file,
        star_type const&)
      :
        out(std::auto_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        pos(0),
        prev_was_string(false),
        io_mode(io_list_directed),
        suppress_new_line_at_end(false)
      {}

      write_loop(
        str_ref const& internal_file,
        char const* fmt)
      :
        out(std::auto_ptr<utils::simple_ostream>(new
          utils::simple_ostream_to_char_ptr_and_size(
            internal_file.elems(), internal_file.len()))),
        internal_file_len(internal_file.len()),
        pos(0),
        prev_was_string(false),
        io_mode(io_formatted),
        fmt_loop(fmt),
        suppress_new_line_at_end(false)
      {}

      std::string const&
      next_edit_descriptor(
        bool final=false)
      {
        while (true) {
          utils::token const* t = fmt_loop.next_executable_token(final);
          if (t == 0) {
            static const std::string empty("");
            return empty;
          }
          std::string const& tv = t->value;
          if (t->type == "string") {
            to_stream_fmt(tv.data(), tv.size());
          }
          else if (t->type == "op") {
            if (tv[0] == ':') {
              // ignored
            }
            else if (tv[0] == '/') {
              to_stream_fmt("\n", 1);
            }
            else if (tv[0] == '$') {
              suppress_new_line_at_end = true;
            }
            else {
              throw BOOST_ADAPTBX_UNREACHABLE_ERROR();
            }
          }
          else if (t->type == "format") {
            if (utils::ends_with_char(tv, 'x')) {
              unsigned n = tv.size();
              if (n != 1) n = utils::unsigned_integer_value(
                tv.data(), n-1);
              for(unsigned i=0;i<n;i++) to_stream(" ", 1);
            }
            else if (std::strchr("adefgilz", tv[0]) != 0) {
              return tv;
            }
            else if (utils::ends_with_char(tv, 'p')) {
              throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 't') {
              throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 's') {
              throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
            }
            else if (tv[0] == 'b') {
              throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
            }
            else {
              throw BOOST_ADAPTBX_UNREACHABLE_ERROR();
            }
          }
          else {
            throw BOOST_ADAPTBX_UNREACHABLE_ERROR();
          }
        }
      }

      write_loop&
      operator,(
        char const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(&val, 1);
        }
        else if (io_mode == io_list_directed) {
          to_stream(&val, 1, /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        char const* val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(val, std::strlen(val));
        }
        else if (io_mode == io_list_directed) {
          to_stream_star(val, std::strlen(val), /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          if (ed[0] == 'a') {
            int n = ed.size();
            ASSERTBX(n+2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt+1, ed.data()+1, n-1);
            fmt[n] = 's';
            fmt[n+1] = '\0';
            char buf[64];
            n = std::sprintf(buf, fmt, val);
            to_stream_fmt(buf, n);
          }
          else {
            to_stream_fmt(val, std::strlen(val));
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        bool const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(bool));
        }
        else if (io_mode == io_list_directed) {
          to_stream_star((val ? "T" : "F"), 1);
          prev_was_string = false;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int n = ed.size();
          if (ed[0] == 'l' && n > 1) {
            n = utils::unsigned_integer_value(ed.data()+1, n-1);
          }
          else {
            n = 1;
          }
          for(int i=1;i<n;i++) to_stream_fmt(" ", 1);
          to_stream_fmt((val ? "T" : "F"), 1);
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_2 const& val)
      {
        if (io_mode == io_unformatted) {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        else if (io_mode == io_list_directed) {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        else {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_4 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_4));
        }
        else if (io_mode == io_list_directed) {
          char buf[64];
          int n = std::sprintf(buf, "%11d", val);
          to_stream_star(buf, n);
          prev_was_string = false;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          if (ed[0] == 'i') {
            int n = ed.size();
            ASSERTBX(n+2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt+1, ed.data()+1, n-1);
            fmt[n] = 'd';
            fmt[n+1] = '\0';
            char buf[64];
            n = std::sprintf(buf, fmt, val);
            to_stream_fmt(buf, n);
          }
          else {
            char buf[64];
            int n = std::sprintf(buf, " %d", val);
            to_stream_fmt(buf, n);
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        integer_star_8 const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(integer_star_8));
        }
        else if (io_mode == io_list_directed) {
          // TODO faster implementation
          std::ostringstream o;
          o.width(21);
          o << val;
          std::string s = o.str();
          to_stream_star(s.data(), s.size());
          prev_was_string = false;
        }
        else {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        float const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(float));
        }
        else if (io_mode == io_list_directed) {
          utils::float_as_string_list_directed conv(val);
          to_stream(conv.begin, conv.n);
          prev_was_string = false;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          if (ed[0] == 'f') {
            int n = ed.size();
            ASSERTBX(n+2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt+1, ed.data()+1, n-1);
            fmt[n] = 'f';
            fmt[n+1] = '\0';
            char buf[64];
            n = std::sprintf(buf, fmt, val);
            to_stream_fmt(buf, n);
          }
          else {
            char buf[64];
            int n = std::sprintf(buf, " %.6g", val);
            to_stream_fmt(buf, n);
          }
        }
        return *this;
      }

      write_loop&
      operator,(
        double const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(
            reinterpret_cast<char const*>(&val),
            sizeof(double));
        }
        else if (io_mode == io_list_directed) {
          utils::double_as_string_list_directed conv(val);
          to_stream(conv.begin, conv.n);
          prev_was_string = false;
        }
        else {
          throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
        }
        return *this;
      }

      write_loop&
      operator,(
        str_cref const& val)
      {
        if (io_mode == io_unformatted) {
          to_stream_unformatted(val.elems(), val.len());
        }
        else if (io_mode == io_list_directed) {
          to_stream(val.elems(), val.len(), /*space*/ !prev_was_string);
          prev_was_string = true;
        }
        else {
          std::string const& ed = next_edit_descriptor();
          int n = ed.size();
          if (ed[0] == 'a' && n > 1) {
            n = utils::unsigned_integer_value(ed.data()+1, n-1);
            to_stream(val.elems(), std::min(val.len(), n));
            for(int i=val.len();i<n;i++) {
              to_stream(" ", 1);
            }
          }
          else {
            to_stream(val.elems(), val.len());
          }
        }
        return *this;
      }

      template <typename T, size_t Ndims>
      write_loop&
      operator,(
        arr_cref<T, Ndims> const& val)
      {
        size_t n = val.size_1d();
        T const* val_begin = val.begin();
        for(size_t i=0;i<n;i++) {
          (*this), val_begin[i];
        }
        return *this;
      }

      template <size_t Ndims>
      write_loop&
      operator,(
        str_arr_cref<Ndims> const& val)
      {
        size_t n = val.size_1d();
        int l = val.len();
        char const* val_begin = val.begin();
        for(size_t i=0;i<n;i++) {
          (*this), str_cref(&val_begin[i*l], l);
        }
        return *this;
      }

      ~write_loop()
      {
        if (out.get() == 0) return;
        if (internal_file_len < 0) {
          if (io_mode == io_unformatted) {
            out->put(end_of_unformatted_record);
            out->put('\0');
          }
          else {
            if (io_mode == io_list_directed) {
              if (pos == 0) out->put(' ');
            }
            else {
              next_edit_descriptor(/*final*/ true);
            }
            if (!suppress_new_line_at_end) {
              out->put('\n');
            }
          }
          out->flush();
        }
        else {
          if (io_mode == io_unformatted) {
            throw BOOST_ADAPTBX_NOT_IMPLEMENTED();
          }
          else {
            if (io_mode == io_list_directed) {
              if (pos == 0) {
                out->put(' ');
                pos++;
              }
            }
            else {
              next_edit_descriptor(/*final*/ true);
            }
            while (pos < internal_file_len) {
              out->put(' ');
              pos++;
            }
          }
        }
      }

    private:

      void
      to_stream(
        char const* buf,
        unsigned n,
        bool space=true)
      {
        switch (io_mode) {
          case io_unformatted:   to_stream_unformatted(buf, n); break;
          case io_list_directed: to_stream_star(buf, n, space); break;
          default:               to_stream_fmt(buf, n);
        }
      }

      void
      to_stream_unformatted(
        char const* buf,
        unsigned n)
      {
        for(unsigned i=0;i<n;i++) {
          char c = buf[i];
          out->put(c);
          if (c == end_of_unformatted_record) {
            out->put(c);
          }
        }
      }

      void
      to_stream_fmt(
        char const* buf,
        unsigned n)
      {
        out->put(buf, n);
      }

      void
      to_stream_star(
        char const* buf,
        unsigned n,
        bool space=true)
      {
        if (pos == 0) {
          out->put(' ');
          pos = 1;
        }
        else if (pos + (space ? 1 : 0) + n > 80) {
          out->put("\n ", 2);
          pos = 1;
        }
        else if (space) {
          out->put(' ');
          pos++;
        }
        out->put(buf, n);
        pos += n;
      }
  };

  struct common_write
  {
    common& cmn;

    common_write(
      common& cmn_)
    :
      cmn(cmn_)
    {}

    write_loop
    operator()(
      int unit,
      unformatted_type const&)
    {
      write_loop result(cmn, unit, unformatted);
      return result;
    }

    write_loop
    operator()(
      int unit,
      star_type const&)
    {
      write_loop result(cmn, unit, star);
      return result;
    }

    write_loop
    operator()(
      int const& unit,
      char const* fmt)
    {
      write_loop result(cmn, unit, fmt);
      return result;
    }

    write_loop
    operator()(
      str_ref const& internal_file,
      star_type const&)
    {
      write_loop result(internal_file, star);
      return result;
    }

    write_loop
    operator()(
      str_ref const& internal_file,
      char const* fmt)
    {
      write_loop result(internal_file, fmt);
      return result;
    }
  };

} // namespace fem

#endif // GUARD
