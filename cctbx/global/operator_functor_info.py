arithmetic_unary_ops = ("-")
arithmetic_binary_ops = ("+", "-", "*", "/", "%")
arithmetic_in_place_binary_ops = ("+=", "-=", "*=", "/=", "%=")

logical_unary_ops = ("!")
logical_binary_ops = ("&&", "||")

boolean_ops = ("==", "!=", ">", "<", ">=", "<=")

std_abinop_function_objects = {
  "+": "std::plus",
  "-": "std::minus",
  "*": "std::multiplies",
  "/": "std::divides",
  "%": "std::modulus",
}

aipbinop_function_objects = {
  "+=": "ip_plus",
  "-=": "ip_minus",
  "*=": "ip_multiplies",
  "/=": "ip_divides",
  "%=": "ip_modulus",
}

std_lbinop_function_objects = {
  "&&": "std::logical_and",
  "||": "std::logical_or",
}

std_boolop_function_objects = {
  "==": "std::equal_to",
  "!=": "std::not_equal_to",
  ">":  "std::greater",
  "<":  "std::less",
  ">=": "std::greater_equal",
  "<=": "std::less_equal",
}

std_all_function_objects = {}
std_all_function_objects.update(std_abinop_function_objects)
std_all_function_objects.update(std_lbinop_function_objects)
std_all_function_objects.update(std_boolop_function_objects)
