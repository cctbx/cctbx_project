      program prog
      implicit integer (a-z)
      arr = 1
      arr_index = arr + 2
      arr_cref = arr_index + 3
      arr_ref = arr_cref + 4
      arr_size = arr_ref + 5
      arr_1d = arr_size + 6
      arr_2d = arr_1d + 7
      arr_3d = arr_2d + 8
      common_read = arr_3d + 9
      common_variant = common_read + 10
      common_write = common_variant + 11
      datum = common_write + 12
      dimension = datum + 13
      dim1 = dimension + 14
      equivalence = dim1 + 15
      local_equivalences = equivalence + 16
      read_loop = local_equivalences + 17
      save_equivalences = read_loop + 18
      star = save_equivalences + 19
      str_arr_cref = star + 20
      str_arr_ref = str_arr_cref + 21
      str_cref = str_arr_ref + 22
      str_index = str_cref + 23
      str_ref = str_index + 24
      values = str_ref + 25
      write_loop = values + 26
      and = write_loop + 27
      and_eq = and + 28
      asm = and_eq + 29
      auto = asm + 30
      bitand = auto + 31
      bitor = bitand + 32
      bool = bitor + 33
      break = bool + 34
      case = break + 35
      catch = case + 36
      char = catch + 37
      class = char + 38
      compl = class + 39
      const = compl + 40
      const_cast = const + 41
      continue = const_cast + 42
      default = continue + 43
      delete = default + 44
      do = delete + 45
      double = do + 46
      dynamic_cast = double + 47
      else = dynamic_cast + 48
      enum = else + 49
      explicit = enum + 50
      export = explicit + 51
      extern = export + 52
      false = extern + 53
      float = false + 54
      for = float + 55
      friend = for + 56
      goto = friend + 57
      if = goto + 58
      inline = if + 59
      int = inline + 60
      long = int + 61
      mutable = long + 62
      namespace = mutable + 63
      new = namespace + 64
      not = new + 65
      not_eq = not + 66
      operator = not_eq + 67
      or = operator + 68
      or_eq = or + 69
      private = or_eq + 70
      protected = private + 71
      public = protected + 72
      register = public + 73
      reinterpret_cast = register + 74
      return = reinterpret_cast + 75
      short = return + 76
      signed = short + 77
      sizeof = signed + 78
      static = sizeof + 79
      static_cast = static + 80
      struct = static_cast + 81
      switch = struct + 82
      template = switch + 83
      this = template + 84
      throw = this + 85
      true = throw + 86
      try = true + 87
      typedef = try + 88
      typeid = typedef + 89
      typename = typeid + 90
      union = typename + 91
      unsigned = union + 92
      using = unsigned + 93
      virtual = using + 94
      void = virtual + 95
      volatile = void + 96
      wchar_t = volatile + 97
      while = wchar_t + 98
      xor = while + 99
      xor_eq = xor + 100
      write(6, *) xor_eq
      call exercise_common_member
      end

      subroutine exercise_common_member
      integer template
      common /vars/ template
      template = 123
      write(6, *) template
      end
