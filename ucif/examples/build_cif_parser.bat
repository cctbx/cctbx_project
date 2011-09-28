cl /DNOMINMAX /EHsc /DANTLR3_NODEBUGGER /Fe./cif_parser ^
/I..\..\ /I ..\..\..\antlr3\include\ /I..\...\..\antlr3 ^
main.cpp ..\cifLexer.cpp ..\cifParser.cpp ..\cifWalker.cpp ^
  "..\..\..\antlr3\src\antlr3baserecognizer.c" ^
  "..\..\..\antlr3\src\antlr3basetree.c" ^
  "..\..\..\antlr3\src\antlr3basetreeadaptor.c" ^
  "..\..\..\antlr3\src\antlr3bitset.c" ^
  "..\..\..\antlr3\src\antlr3collections.c" ^
  "..\..\..\antlr3\src\antlr3commontoken.c" ^
  "..\..\..\antlr3\src\antlr3commontree.c" ^
  "..\..\..\antlr3\src\antlr3commontreeadaptor.c" ^
  "..\..\..\antlr3\src\antlr3commontreenodestream.c" ^
  "..\..\..\antlr3\src\antlr3convertutf.c" ^
  "..\..\..\antlr3\src\antlr3cyclicdfa.c" ^
  "..\..\..\antlr3\src\antlr3debughandlers.c" ^
  "..\..\..\antlr3\src\antlr3encodings.c" ^
  "..\..\..\antlr3\src\antlr3exception.c" ^
  "..\..\..\antlr3\src\antlr3filestream.c" ^
  "..\..\..\antlr3\src\antlr3inputstream.c" ^
  "..\..\..\antlr3\src\antlr3intstream.c" ^
  "..\..\..\antlr3\src\antlr3lexer.c" ^
  "..\..\..\antlr3\src\antlr3parser.c" ^
  "..\..\..\antlr3\src\antlr3rewritestreams.c" ^
  "..\..\..\antlr3\src\antlr3string.c" ^
  "..\..\..\antlr3\src\antlr3stringstream.c" ^
  "..\..\..\antlr3\src\antlr3tokenstream.c" ^
  "..\..\..\antlr3\src\antlr3treeparser.c" ^
  "..\..\..\antlr3\src\antlr3ucs2inputstream.c"
