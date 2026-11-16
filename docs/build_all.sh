#!/usr/bin/env bash
# Rebuild figures/tables and all LaTeX documents (ES/EN).
set -e
cd "$(dirname "$0")/.."
python3 examples/generate_all.py
cd docs
for d in manual:manual_es manual:manual_en article:article_es article:article_en beamer:beamer_es beamer:beamer_en course:curso_es course:course_en; do
  dir=${d%%:*}; f=${d##*:}
  (cd "$dir" && latexmk -pdf -interaction=nonstopmode -halt-on-error "$f.tex" > /dev/null && echo "built $dir/$f.pdf")
done
