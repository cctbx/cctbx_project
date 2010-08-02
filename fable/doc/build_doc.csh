#! /bin/csh -f
if ($#argv != 1) then
  echo "Missing command-line argument: target_dir"
  exit 1
endif
set verbose
set target_dir="$1"
set fable_dist="`libtbx.show_dist_paths fable`"
#
cp "$fable_dist/../dox/rst/default.css" "$target_dir"
#
fable.cout "$fable_dist/test/valid/sf.f" --no-fem-do-safe --namespace=example > "$target_dir/sf.cpp"
#
fable.cout "$fable_dist/test/valid/conv_recipe.f" --no-fem-do-safe --namespace=conv_recipe > "$target_dir/conv_recipe.cpp"
#
fable.cout "$fable_dist/test/valid/dp_example.f" --no-fem-do-safe --namespace=dp_example --dynamic-parameter="int dp_example=100" > "$target_dir/dp_example.cpp"
#
fable.cout "$fable_dist/test/valid/common_variants.f" --no-fem-do-safe --namespace=common_variants > "$target_dir/common_variants.cpp"
mv fable_cout_common_report "$target_dir"
#
(cd "$fable_dist/.."; tar cf "$target_dir/tmp.tar" --exclude='.svn' --exclude='*.pyc' libtbx tbxx fable)
(cd "$target_dir" && rm -rf sources; mkdir sources; cd sources; tar xf ../tmp.tar; rm ../tmp.tar; mv tbxx fable; cp -a "`libtbx.find_in_repositories lapack_fem`" .)
#
cp index.txt "$target_dir"
docutils.rst2html -stg --report=error index.txt > "$target_dir/index.html"
