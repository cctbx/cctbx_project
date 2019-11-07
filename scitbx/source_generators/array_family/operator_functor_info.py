from __future__ import absolute_import, division, print_function
arithmetic_unary_ops = ("-")
arithmetic_binary_ops = ("+", "-", "*", "/", "%")
arithmetic_in_place_binary_ops = ("+=", "-=", "*=", "/=", "%=")

logical_unary_ops = ("!")
logical_binary_ops = ("&&", "||")

boolean_binary_ops = ("==", "!=", ">", "<", ">=", "<=")

unary_functors = {
  "-": "negate",
  "!": "logical_not",
}

binary_functors = {
  "+": "plus",
  "-": "minus",
  "*": "multiplies",
  "/": "divides",
  "%": "modulus",
  "&&": "logical_and",
  "||": "logical_or",
  "==": "equal_to",
  "!=": "not_equal_to",
  ">":  "greater",
  "<":  "less",
  ">=": "greater_equal",
  "<=": "less_equal",
}

in_place_binary_functors = {
  "+=": "ip_plus",
  "-=": "ip_minus",
  "*=": "ip_multiplies",
  "/=": "ip_divides",
  "%=": "ip_modulus",
}
