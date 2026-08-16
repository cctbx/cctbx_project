#!/bin/sh
# Change the version label across the updated docs, if v121 is not right.
#
#     sh RENUMBER.sh <docs_dir> <new_label>
#     sh RENUMBER.sh $LANGCHAIN_DIR/docs v120.5
#
# 27 occurrences across 4 files. The label was chosen by me, not by Tom:
# the last released is v120.4, and whether this is 121 or 120.5 is a
# judgement about scope that I had no basis to make.
set -e
D="$1"; NEW="$2"
if [ -z "$D" ] || [ -z "$NEW" ]; then
  echo "usage: sh RENUMBER.sh <docs_dir> <new_label>"; exit 1
fi
n=0
for f in "$D"/ARCHITECTURE.md "$D"/OVERVIEW.md "$D"/CHANGELOG.md \
         "$D"/AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md ; do
  [ -f "$f" ] || continue
  c=`grep -o 'v121' "$f" | wc -l | tr -d ' '`
  if [ "$c" != "0" ] ; then
    sed -i.bak "s/v121/$NEW/g" "$f" && rm -f "$f.bak"
    echo "  $f: $c replaced"
    n=`expr $n + $c`
  fi
done
echo "  total $n occurrences -> $NEW"
echo "  also check: '## Version 121 (' heading and the CHANGELOG title line"
