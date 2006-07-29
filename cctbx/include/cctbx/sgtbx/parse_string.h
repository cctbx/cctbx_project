#ifndef CCTBX_SGTBX_PARSE_STRING_H
#define CCTBX_SGTBX_PARSE_STRING_H

#include <cstddef>
#include <string>

namespace cctbx { namespace sgtbx {

  //! Class for communicating string parsing errors.
  /*! This class is used by functions such as
      cctbx::sgtbx::space_group::parse_hall_symbol()
      or a constructor of class cctbx::sgtbx::rt_mx to communitcate errors
      that are detected during the interpretation of an input string.<br>
      Intended use:<pre>
      using namespace cctbx::sgtbx;
      parse_string hall_symbol("P x");
      space_group sg;
      try {
        sg.parse_hall_symbol(hall_symbol);
      }
      catch (cctbx::error const& e) {
        std::cout << e.what() << std::endl;
        std::cout << hall_symbol.string() << std::endl;
        for (std::size_t i = 0; i < hall_symbol.where(); i++) std::cout << "_";
        std::cout << "^" << std::endl;
      }
      </pre>
      This will produce the output:<pre>
      cctbx Error: Improper symbol for rotational order.
      P x
      __^
      </pre>
   */
  class parse_string
  {
    public:
      //! Initializes the parse_string with the empty string.
      parse_string()
      :
        pos_(0), marked_pos_(0)
      {}

      //! Initializes the parse_string with the input string.
      explicit
      parse_string(std::string const& str)
      :
        s_(str
#if defined(__GNUC__) && __GNUC__ == 3 \
                      && (__GNUC_MINOR__ == 2 || __GNUC_MINOR__ == 3)
              +"" // work around gcc 3.2 & 3.3 bug
#endif
                  ),
        pos_(0),
        marked_pos_(0)
      {}

      //! Returns the input string.
      const char*
      string() const { return s_.c_str(); }

      //! Index of last character accessed by the parsing algorithm.
      std::size_t
      where() const { return pos_; }

      //! Support for formatting exceptions.
      /*! Produces output of the form:
            P x
            __^
       */
      std::string
      format_where_message(std::string const& prefix) const
      {
        std::string result = prefix + s_ + "\n" + prefix;
        for(std::size_t i=0;i<pos_;i++) result += "_";
        result += "^";
        return result;
      }

      //! For internal use only.
      char
      operator()() const { return s_[pos_]; }

      //! For internal use only.
      const char*
      peek() { return s_.c_str() + pos_; }

      //! For internal use only.
      char
      get()
      {
        if (pos_ >= s_.size()) return '\0';
        return s_[pos_++];
      }

      //! For internal use only.
      void
      skip(std::size_t n = 1) { pos_ += n; }

      //! For internal use only.
      void
      set_mark() { marked_pos_ = pos_; }

      //! For internal use only.
      void
      go_to_mark() { pos_ = marked_pos_; }

    private:
      std::string s_;
      std::size_t pos_;
      std::size_t marked_pos_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_PARSE_STRING_H
