dox=`libtbx.find_in_repositories dox`
TEXINPUTS=".:$dox:" exec pdflatex "$@"
