#! /bin/csh -fe

if (! -d cctbx_project) then
  echo "ERROR: cctbx_project not found in current working directory."
  exit 1
endif
mkdir FROM_BUNDLE
if (-f TAG) then
  echo "mv TAG FROM_BUNDLE"
  mv TAG FROM_BUNDLE
endif
foreach f (*)
  if (-d "$f") then
    if (-d "cctbx_project/$f") then
      echo "mv $f FROM_BUNDLE"
            mv "$f" FROM_BUNDLE
    endif
  endif
end
exit 0
