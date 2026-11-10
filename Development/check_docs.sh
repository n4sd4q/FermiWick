#!/bin/sh
# Checks that every command of the documents gives what the documents show, in fresh kernels:
#  - the showcase and the manual: their examples are run again and every output must equal the one in examples/;
#  - the guided tour Examples.wl runs as a script without messages;
#  - the first example of the README gives the formulas printed after its code;
#  - every page of the wiki (wiki/*.md) gives what it shows after each block (Development/check_wiki.wl), the pages at the
#    same time.
# Run from anywhere:  sh Development/check_docs.sh      (exit status 0 when everything passes)
cd "$(dirname "$0")/.." || exit 1
status=0
tmp=$(mktemp -d)
for doc in Showcase Documentation; do
  SRC=$doc/examples.wl OUT=$tmp/$doc math -script Development/run_examples.wl > $tmp/$doc.log 2>&1
  n=$(ls $doc/examples | wc -l)
  if diff -r -q $doc/examples $tmp/$doc > $tmp/$doc.diff; then echo "ok   $doc: $n files, all as published"
  else echo "FAIL $doc: outputs differ from the published ones:"; sed 's/^/     /' $tmp/$doc.diff; status=1; fi
done
out=$(math -script Examples.wl 2>&1 | grep -v '^$')
if [ -z "$out" ]; then echo "ok   Examples.wl: runs without messages"; else echo "FAIL Examples.wl:"; echo "$out" | head -20; status=1; fi
math -script Development/check_readme.wl | grep -v '^$' | sed 's/^/     /' || status=1
math -script Development/check_readme.wl > /dev/null 2>&1 || status=1
for page in wiki/*.md; do math -script Development/check_wiki.wl "$page" > "$tmp/wiki_$(basename "$page").log" 2>&1 & done; wait
for page in wiki/*.md; do log="$tmp/wiki_$(basename "$page").log"; tail -1 "$log" | sed 's/^/     /'
  grep -q '^ok ' "$log" || { grep -v '^$' "$log" | head -20 | sed 's/^/     /'; status=1; }; done
python3 Development/make_site.py --check | sed 's/^/ok   /' || status=1
python3 Development/make_site.py --check > /dev/null 2>&1 || status=1
rm -rf $tmp
exit $status
