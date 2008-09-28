dox=`libtbx.find_in_repositories dox`
BIBINPUTS=".:$dox" exec bibtex "$@"
