#! /bin/csh -f
set verbose
set target_dir="$HOME/public_html/fable"
set fable_dist="`libtbx.show_dist_paths fable`"
cp "$fable_dist/test/valid/sf.f" "$target_dir"
fable.cout "$fable_dist/test/valid/sf.f" --no-fem-do-safe --namespace=example > "$target_dir/sf.cpp"
cp "`libtbx.find_in_repositories lapack_fem/dsyev.hpp`" "$target_dir"
docutils_dev.rst2html -stg index.txt > "$target_dir/index.html"
