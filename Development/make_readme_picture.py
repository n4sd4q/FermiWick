"""Draws the picture of Wick's theorem of the README, Documentation/readme_wick.svg, from Documentation/readme/wick.tex
(the output of WickTeX computed for the showcase): pdflatex, pdftocairo -svg, then the ink in currentColor on a
transparent page, with the colour of GitHub's text in its light and dark themes, so that the picture reads as the
formulas around it. The website inlines the same file without that style, and its ink takes the colour of the page.
    python3 Development/make_readme_picture.py"""
import os, re, subprocess, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "Documentation", "readme")
OUT = os.path.join(ROOT, "Documentation", "readme_wick.svg")
STYLE = ('<style>svg { color: #1f2328; } @media (prefers-color-scheme: dark) { svg { color: #e6edf3; } }</style>')

with tempfile.TemporaryDirectory() as tmp:
    for run in range(2):   # simpler-wick places the contractions on the second run (TikZ, remember picture)
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", tmp, "wick.tex"], cwd=SRC, check=True,
                       capture_output=True)
    subprocess.run(["pdftocairo", "-svg", os.path.join(tmp, "wick.pdf"), os.path.join(tmp, "wick.svg")], check=True)
    svg = open(os.path.join(tmp, "wick.svg"), encoding="utf-8").read()
svg = re.sub(r'(fill|stroke)="rgb\(0%, 0%, 0%\)"', r'\1="currentColor"', svg)
if re.search(r'rgb\(\d', svg):
    raise SystemExit("a colour other than black is left in the picture")
svg = svg.replace(">", ">\n" + STYLE, 1) if svg.startswith("<svg") else re.sub(r"(<svg[^>]*>)", r"\1\n" + STYLE, svg, count=1)
open(OUT, "w", encoding="utf-8").write(svg)
print("wrote", os.path.relpath(OUT, ROOT))
