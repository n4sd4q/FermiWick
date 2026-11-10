#!/bin/sh
# Computes the examples and typesets the showcase.
cd "$(dirname "$0")"
export TEXINPUTS="../Development:$TEXINPUTS"   # the common style, n4sd4q.sty
SRC=examples.wl OUT=examples math -script ../Development/run_examples.wl > examples.log 2>&1 || exit 1
pdflatex -interaction=nonstopmode -halt-on-error FermiWick_Showcase.tex > /dev/null && pdflatex -interaction=nonstopmode FermiWick_Showcase.tex > /dev/null
grep -E "Warning|Overfull|Error|missing example" FermiWick_Showcase.log | grep -v "Font Warning" | head
rm -f FermiWick_Showcase.aux FermiWick_Showcase.out FermiWick_Showcase.toc
