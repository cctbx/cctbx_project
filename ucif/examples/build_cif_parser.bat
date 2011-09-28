SET LIBTBX_ANTLR3=..\..\..\antlr3
cl /DNOMINMAX /EHsc /DANTLR3_NODEBUGGER /Fe./cif_parser ^
/I..\..\ /I %LIBTBX_ANTLR3%\include\ /I%LIBTBX_ANTLR3%\ ^
main.cpp ..\cifLexer.cpp ..\cifParser.cpp ..\cifWalker.cpp ^
  "%LIBTBX_ANTLR3%\src\antlr3baserecognizer.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3basetree.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3basetreeadaptor.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3bitset.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3collections.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3commontoken.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3commontree.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3commontreeadaptor.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3commontreenodestream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3convertutf.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3cyclicdfa.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3debughandlers.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3encodings.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3exception.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3filestream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3inputstream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3intstream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3lexer.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3parser.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3rewritestreams.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3string.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3stringstream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3tokenstream.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3treeparser.c" ^
  "%LIBTBX_ANTLR3%\src\antlr3ucs2inputstream.c"
