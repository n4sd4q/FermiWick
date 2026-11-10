#!/usr/bin/env python3
"""The demo video of the package: Media/demo.mp4 (1920 x 1080, 60 fps, H.264; its first frame, the one a player shows
before it plays, is the title card in full view), its looping preview for the README,
Media/demo_preview.webp (960 x 540, 60 fps; its first frame, the one shown where animated images do not play, is a card
in full view), and its poster, Media/demo_poster.png (that title card: the website shows it until the video plays).

Nothing in it is written by hand: the session (CARDS) is made of the code blocks of "A first example" of README.md, each
with what the README prints after it (the formulas, which Development/check_docs.sh keeps equal to what ToTeX writes,
with GitHub's escapes undone; the picture of Wick's theorem is Documentation/readme/wick.tex, the output of WickTeX
typeset with simpler-wick), and of examples computed for the showcase and the manual (Showcase/examples and
Documentation/examples: the code NAME.wl, the output NAME.tex, NAME.txt or, for WickTeX, NAME_typeset.tex). One card for
each thing the package does, from the ladder operators to LaTeX: a heading, a second line naming the examples, and a few
cells typed one after the other as in a notebook (In[n]:= and Out[n]= numbered as one session numbers them, the cell
brackets on the right, Mathematica's colours), each followed by its output. The cells of CARDS, in their order, are one
session: run in one fresh kernel, each gives the output shown (the variables a cell uses are set by a cell before it;
--check runs it). The pictures go through pdflatex (twice: simpler-wick places its contractions on the second run) and
pdftocairo (one PNG per element, transparent), the frames are put together with PIL (at 1920 x 1080; the preview is
reduced from them) and encoded by ffmpeg. The logo of
the first and the last card is drawn again in every frame from the geometry of Logo/FermiWick_logo.svg (read from the
file): on the title card its arrow heads travel once along their contractions, from the creator (the open circle) into
the annihilator (the filled circle) and out of the creator again, with a soft start, a glide and a slight overshoot
before they settle; the closing card is the title card again (its last frame is the first frame of the video, so that
the loop has no seam), its logo drawing itself: the four operators appear, then each contraction runs from its creator
to its annihilator, its arrow head riding behind the pen and settling on its place. Needs pdflatex (TeX Live:
mathpazo, Source Code Pro, simpler-wick), pdftocairo (poppler), python3 with PIL, ffmpeg with libx264. Run from any
folder:
    python3 Development/make_demo_video.py            (all three files, about eight minutes: most of it the WebP)
    python3 Development/make_demo_video.py --frames   (also a few frames as PNG in build/video/frames, to look at;
                                                      with --only-frames, nothing else)
    python3 Development/make_demo_video.py --poster   (only the poster, a few seconds)
    python3 Development/make_demo_video.py --check    (only the check of the session, in Mathematica: wolfram -script,
                                                      a minute or so)
"""
import concurrent.futures
import functools
import math
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(ROOT, "build", "video")
MEDIA = os.path.join(ROOT, "Media")
W, H, FPS = 1920, 1080, 60
# Every length in pixels below (fonts, margins, gaps, line widths, brackets) is written for a frame 1280 px wide and
# scaled by S, so that a frame of another size is the same picture, only sharper
S = W / 1280


def px(v):
    """A length in pixels of a frame 1280 px wide, in pixels of this frame."""
    return v * S


def ipx(v):
    """The same, in whole pixels."""
    return int(math.floor(v * S + 0.5))


# The preview: an animated WebP, which can run at 60 fps (a GIF cannot: its delays are in hundredths of a second, so
# 50 fps at most, and browsers slow down delays below 0.02 s), in full colour; GitHub shows it in the README as a plain
# <img>, and it plays wherever WebP does (every current browser, the GitHub apps)
PREVIEW_W, PREVIEW_H, PREVIEW_FPS, PREVIEW_QUALITY = 960, 540, 60, 90
# the LaTeX pages are drawn for a frame 16 cm wide: 120 px to the centimetre
DPI = W / (16 / 2.54)
PT = DPI / 72.27        # pixels to the TeX point
BORDER = 3.0            # the border of every LaTeX element (pt, the option of standalone)

# the colours of the package ("Orchid", Development/n4sd4q-palette.tex): the solid accent of the titles and the band;
# those of the logo (the gradient of its badge, white and the light tint) are read from Logo/FermiWick_logo.svg
ACCENT = (138, 40, 128)
MUTED = (89, 99, 110)
WHITE = (255, 255, 255)
# Mathematica's colours in an input cell: comments and strings grey, undefined symbols blue, the variables of a Sum
# green, the rest (the system's symbols, the package's, those already given a value) black
CODE_COLORS = {"comment": (140, 140, 140), "string": (110, 110, 110), "symbol": (36, 74, 170), "local": (52, 128, 84),
               "plain": (0, 0, 0)}
LABEL = (70, 96, 140)        # In[n]:= and Out[n]=
BRACKET = (150, 172, 200)    # the cell brackets on the right


# ---------------------------------------------------------------------------------------------------------------------
# the example of the README

def readme_section():
    text = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    parts = re.split(r"^## A first example\s*$", text, maxsplit=1, flags=re.M)
    if len(parts) < 2:
        sys.exit('README.md has no section "## A first example"')
    return re.split(r"^## ", parts[1], maxsplit=1, flags=re.M)[0]


def from_github(tex):
    """A formula of the README as LaTeX: the escapes it needs on GitHub (whose Markdown runs before MathJax) undone."""
    tex = tex.replace("\\\\{", "\\{").replace("\\\\}", "\\}").replace("^{\\ast}", "^*")
    tex = re.sub(r"\\thinspace(?![A-Za-z])", r"\\,", tex)
    tex = re.sub(r"\\lt(?![A-Za-z])", "<", tex)
    return re.sub(r"\\gt(?![A-Za-z])", ">", tex).strip()


def wick_picture():
    """The formula of Documentation/readme/wick.tex, the picture of Wick's theorem in the README, its \\input read in."""
    path = os.path.join(ROOT, "Documentation", "readme", "wick.tex")
    text = open(path, encoding="utf-8").read()
    body = text.split("\\begin{document}", 1)[1].split("\\end{document}", 1)[0].strip()

    def read(m):
        name = os.path.join(os.path.dirname(path), m.group(1))
        name = name if os.path.exists(name) else name + ".tex"
        return open(name, encoding="utf-8").read().strip()
    return re.sub(r"\\input\{([^}]+)\}", read, body)


def readme_steps():
    """The code blocks of "A first example", each with what the README prints after it: [code, LaTeX or None]."""
    steps = []
    for m in re.finditer(r"```mathematica\n(.*?)```|^\$\$([^\n]*?)\$\$[ \t]*$|<img src=\"([^\"]+)\"", readme_section(),
                         re.S | re.M):
        if m.group(1) is not None:
            steps.append([m.group(1).strip("\n"), None])
        elif steps and steps[-1][1] is None:
            if m.group(2) is not None:
                steps[-1][1] = from_github(m.group(2))
            elif m.group(3).endswith("readme_wick.svg"):
                steps[-1][1] = wick_picture()
    # the loading line, then the seven steps, each with its output
    if len(steps) != 8 or steps[0][1] is not None or any(out is None for _, out in steps[1:]):
        sys.exit('the section "A first example" of README.md is not the one this script knows: the loading line, then '
                 "seven code blocks, each followed by its formula (the last by readme_wick.svg)")
    return steps


EXAMPLES = {"showcase": os.path.join(ROOT, "Showcase", "examples"), "manual": os.path.join(ROOT, "Documentation", "examples")}


def example(doc, name):
    """An example computed for the showcase or the manual (Development/run_examples.wl): its code, NAME.wl, and its output,
    ("tex", the LaTeX of NAME.tex, or of NAME_typeset.tex for a string of LaTeX that is a formula: WickTeX) or ("text",
    NAME.txt), or None."""
    path = os.path.join(EXAMPLES[doc], name)
    if not os.path.exists(path + ".wl"):
        sys.exit("missing %s.wl: run the examples of the %s first" % (path, doc))
    code = open(path + ".wl", encoding="utf-8").read().strip("\n")
    for ext, kind in (("_typeset.tex", "tex"), (".tex", "tex"), (".txt", "text")):
        if os.path.exists(path + ext):
            return code, (kind, open(path + ext, encoding="utf-8").read().strip())
    return code, None


# The session: one card for each thing the package does, its heading and a second line naming the examples, and its
# cells, in the order of one session (the variables a cell uses are set by a cell before it). Nothing is typed here: a
# cell is ("readme", k), the code block k of "A first example" with what the README prints after it (0 is the line that
# loads the package), or ("showcase", NAME) or ("manual", NAME), an example computed for the showcase or the manual,
# code and output (example). The last number is how long the card stays in full view (s).
CARDS = [
    ("Ladder operators and normal order", "anticommutation, the Pauli principle, sites and spins as symbols",
     [("readme", 0), ("readme", 1), ("showcase", "pauli"), ("showcase", "samesite"), ("showcase", "spinless")], 2.0),
    ("Spins from ladder operators", "the exchange identity, spin flips, the chirality from a commutator",
     [("readme", 2), ("showcase", "decompose"), ("showcase", "spinlength"), ("showcase", "commutator")], 2.0),
    ("Exact diagonalisation", "any model as on paper: the Hubbard dimer, a ring of eight sites",
     [("readme", 3), ("showcase", "ring8")], 2.6),
    ("Wick's theorem in the state you assume", "the energy; the chirality with and without time reversal",
     [("readme", 4), ("readme", 5), ("manual", "symmetries")], 2.4),
    ("Mean field of any order", "the Hubbard term in a magnetic state, the chirality with pairing",
     [("readme", 6), ("showcase", "singletstate"), ("showcase", "chiralfields")], 2.4),
    ("Gutzwiller projectors", "the projector, the double occupancy, the hopping renormalised",
     [("manual", "load"), ("showcase", "gutzwillerprojector"), ("showcase", "gutzwillerd"), ("showcase", "hopping")], 2.2),
    ("LaTeX", "any result as LaTeX; Wick's theorem written out and drawn",
     [("showcase", "totex"), ("readme", 7), ("manual", "wicktex")], 2.8),
]


SESSION_CHECK = r'''(* written by Development/make_demo_video.py --check: the cells of the video, in their order, in one fresh
   kernel; each result is compared with the output the video shows, written as the examples write it (ToTeX with the
   replacements of Development/run_examples.wl, the text in InputForm) *)
PacletDirectoryLoad[FileNameJoin[{"%(root)s", "FermiWick"}]];
docPost[t_String] := StringReplace[t, {
    RegularExpression["\\\\text\\{Hop\\}\\(([^()]*)\\)"] :> "T_{" <> StringReplace[StringDelete["$1", " "], "," -> ""] <> "}",
    RegularExpression["(?<![A-Za-z\\\\])m\\(([^()]*)\\)"] :> "m_{$1}",
    RegularExpression["\\\\text\\{nk\\}\\(([^()]*)\\)"] :> "n_{$1}",
    "\\text{eps0}" -> "\\varepsilon_0", "\\text{Jc}" -> "J_\\chi"}];
toText[s_String] := s;
toText[e_] := StringReplace[ToString[e /. FermiWick`GaussianState[a_Association] :> stateShown[Lookup[a, "Description", "..."]],
    InputForm, PageWidth -> 110], "stateShown" -> "GaussianState"];
read[name_] := Import[FileNameJoin[{"%(dir)s", name}], "Text", CharacterEncoding -> "UTF8"];
bad = 0;
Do[
  {name, kind, source} = cell;
  result = ToExpression[read[name <> ".wl"]];
  If[kind =!= "none",
   got = StringTrim[Switch[kind, "tex", docPost[FermiWick`ToTeX[result]], "string", docPost[result], "text", toText[result]]];
   If[got === StringTrim[read[name <> ".out"]], Print["ok    ", source],
    bad++; Print["DIFF  ", source, "\n  computed: ", got, "\n  shown:    ", read[name <> ".out"]]]],
  {cell, %(cells)s}];
Print[If[bad == 0, "every output as the video shows it", ToString[bad] <> " outputs differ"]];
Exit[If[bad == 0, 0, 1]]
'''


def check_session():
    """--check: the session of CARDS run in one fresh kernel (wolfram -script), each result compared with the output the
    video shows (SESSION_CHECK, written to build/video/session with one file for each cell and its output)."""
    folder = os.path.join(BUILD, "session")
    os.makedirs(folder, exist_ok=True)
    cells, readme = [], readme_steps()
    for k, (_, _, refs, _) in enumerate(CARDS, 1):
        for j, (doc, key) in enumerate(refs):
            code, out = readme[key] if doc == "readme" else example(doc, key)
            if doc == "readme" and out:
                # the README's formulas are ToTeX of the result; its picture, the typeset WickTeX of the showcase
                kind = "string" if out.startswith("$") else "tex"
                out = example("showcase", "wicktex")[1][1] if kind == "string" else out
            elif out:
                kind = "string" if os.path.exists(os.path.join(EXAMPLES[doc], key + "_typeset.tex")) else out[0]
                out = out[1]
            else:
                kind = "none"
            name = "cell%02d_%d" % (k, j)
            with open(os.path.join(folder, name + ".wl"), "w", encoding="utf-8") as fh:
                fh.write(code)
            if out:
                with open(os.path.join(folder, name + ".out"), "w", encoding="utf-8") as fh:
                    fh.write(out)
            cells.append('{"%s", "%s", "%s %s"}' % (name, kind, doc, key))
    path = os.path.join(folder, "session_check.wl")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(SESSION_CHECK % dict(root=ROOT, dir=folder, cells="{" + ", ".join(cells) + "}"))
    print("the session in", os.path.relpath(path, ROOT), flush=True)
    if not shutil.which("wolfram"):
        sys.exit("wolfram not found: run  wolfram -script %s" % path)
    sys.exit(subprocess.run(["wolfram", "-script", path]).returncode)


def session_cells():
    """The cards of CARDS with their cells read: [(heading, second line, [(code, output)], hold)]."""
    readme = readme_steps()
    out = []
    for head, line, cells, hold in CARDS:
        read = []
        for doc, key in cells:
            if doc == "readme":
                code, tex = readme[key]
                read.append((code, ("tex", tex) if tex else None))
            else:
                read.append(example(doc, key))
        out.append((head, line, read, hold))
    return out


# ---------------------------------------------------------------------------------------------------------------------
# the pictures: LaTeX elements, one page each

PREAMBLE = r"""\documentclass[12pt,multi=elem,border=%gpt]{standalone}
\usepackage[T1]{fontenc}
\usepackage[sc]{mathpazo}
\usepackage[scale=0.9]{sourcecodepro}
\usepackage{amsmath,xcolor,simpler-wick}
\definecolor{accent}{HTML}{%02X%02X%02X}\definecolor{muted}{HTML}{%02X%02X%02X}
\newsavebox{\elembox}
\begin{document}
""" % ((BORDER,) + ACCENT + MUTED)

# the outputs and the pictures of Wick's theorem, a little larger than the code (CODE_SIZE px), as a notebook shows them
OUTPUT_SIZE = r"\fontsize{8.5}{10.5}\selectfont"
CODE_SIZE = px(18)


def output_name(k, j):
    return "o%d_%d" % (k, j)


def elements(cards):
    """name -> LaTeX of every picture of the video."""
    e = {}
    # title and end card
    e["name"] = r"{\color{accent}\fontsize{30}{36}\selectfont\bfseries FermiWick}"
    e["tagline"] = r"{\fontsize{16}{20}\selectfont Fermion operators, as symbols and as matrices}"
    e["url"] = r"{\color{accent}\fontsize{15}{19}\selectfont\ttfamily n4sd4q.github.io/FermiWick}"
    # the captions of the cards and the outputs in LaTeX, in display style (as $$...$$ in the README)
    for k, (head, line, cells, _) in enumerate(cards, 1):
        e["h%d" % k] = r"{\color{accent}\fontsize{19}{23}\selectfont\bfseries %s}" % head
        e["s%d" % k] = r"{\color{muted}\fontsize{13}{16}\selectfont %s}" % line
        for j, (_, out) in enumerate(cells):
            if out and out[0] == "tex":
                body = out[1] if out[1].startswith("$") else r"$\displaystyle %s$" % out[1]     # wick.tex: a formula already
                e[output_name(k, j)] = r"{%s%s}" % (OUTPUT_SIZE, body)
    return e


def text_output(text, fnt):
    """An output written as text (a string, a list of numbers, a state), as a notebook shows it: in the font of the code,
    on a transparent image with the border of the LaTeX elements; and the height of its first baseline."""
    lines, border = text.split("\n"), BORDER * PT
    ascent, descent = fnt.getmetrics()
    lh = px(26)
    img = Image.new("RGBA", (int(max(fnt.getlength(l) for l in lines) + 2 * border + 1),
                             int(2 * border + ascent + descent + lh * (len(lines) - 1) + 1)), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    for k, l in enumerate(lines):
        d.text((border, border + ascent + k * lh), l, font=fnt, fill=(0, 0, 0, 255), anchor="ls")
    return img, border + ascent


def render_elements(cards):
    """Typeset every element and return (name -> RGBA image, transparent, 120 px to the cm; name -> the height of its
    baseline from the top of the image, px)."""
    os.makedirs(BUILD, exist_ok=True)
    els = elements(cards)
    names = list(els)
    with open(os.path.join(BUILD, "elements.tex"), "w", encoding="utf-8") as fh:
        fh.write(PREAMBLE)
        for name in names:
            # each element in a box, its height and depth written to the log: where its baseline is
            fh.write("%% %s\n\\begin{elem}\\sbox{\\elembox}{%s}\\typeout{ELEMBOX %s \\the\\ht\\elembox\\space\\the\\dp\\elembox}"
                     "\\usebox{\\elembox}\\end{elem}\n" % (name, els[name], name))
        fh.write("\\end{document}\n")
    env = dict(os.environ, max_print_line="100000")
    for _ in range(2):      # simpler-wick draws its contractions where the first run has put the operators
        run = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "elements.tex"], cwd=BUILD,
                             capture_output=True, text=True, env=env)
        if run.returncode:
            sys.exit("pdflatex failed:\n" + run.stdout[-3000:])
    baselines = {m.group(1): (BORDER + float(m.group(2))) * PT
                 for m in re.finditer(r"^ELEMBOX (\S+) ([0-9.]+)pt ([0-9.]+)pt$", run.stdout, re.M)}
    pngs = os.path.join(BUILD, "png")
    os.makedirs(pngs, exist_ok=True)
    images = {}

    def page(i, name):
        subprocess.run(["pdftocairo", "-png", "-transp", "-r", "%.4f" % DPI, "-f", str(i), "-l", str(i), "-singlefile",
                        os.path.join(BUILD, "elements.pdf"), os.path.join(pngs, name)], check=True)

    # one pdftocairo per page, several at a time
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count() or 1) as pool:
        list(pool.map(page, range(1, len(names) + 1), names))
    for name in names:
        images[name] = Image.open(os.path.join(pngs, name) + ".png").convert("RGBA")
        box = images[name].getchannel("A").getbbox()
        if box and (box[0] == 0 or box[1] == 0 or box[2] == images[name].width or box[3] == images[name].height):
            print("warning: element %s is cut at the edge of its page" % name)
    # the outputs written as text, in the font of the code (regular, as a notebook writes its outputs)
    textfont = font("SourceCodePro-Regular.otf", CODE_SIZE)
    for k, (_, _, cells, _) in enumerate(cards, 1):
        for j, (_, out) in enumerate(cells):
            if out and out[0] == "text":
                images[output_name(k, j)], baselines[output_name(k, j)] = text_output(out[1], textfont)
    return images, baselines


# ---------------------------------------------------------------------------------------------------------------------
# code, typed as in a notebook

def font(name, size):
    path = subprocess.run(["kpsewhich", name], capture_output=True, text=True).stdout.strip()
    if not path:
        sys.exit("font %s not found (TeX Live package sourcecodepro)" % name)
    return ImageFont.truetype(path, size)


def package_names():
    """The public names of the package (those with a usage message), black in an input cell once it is loaded."""
    names = set()
    kernel = os.path.join(ROOT, "FermiWick", "Kernel")
    for f in sorted(os.listdir(kernel)):
        if f.endswith(".wl"):
            names |= set(re.findall(r"^([A-Za-z$][A-Za-z0-9$]*)::usage", open(os.path.join(kernel, f), encoding="utf-8").read(), re.M))
    return names


NAME = r"[A-Za-z$\u0370-\u03ff][A-Za-z0-9$\u0370-\u03ff]*"
# the names of the system, black in an input cell (those of the examples of the documents, and a few more; another one
# would be coloured as undefined: add it here)
SYSTEM_NAMES = set("""Abs Apply Block C Chop Conjugate Cos D Det E Eigensystem Eigenvalues Exp Expand Factor False First Flatten
    FullSimplify Function I IdentityMatrix If Join K KroneckerDelta Last Length Map MatrixForm Max Min Mod Module N Needs None
    Normal O Pi Print Range Round Row Select Simplify Sin Solve Sort Sqrt Subscript Sum Table Thread Total Transpose True
    With""".split())
TOKEN = re.compile(r"\(\*.*?\*\)|\"(?:[^\"\\]|\\.)*\"|" + NAME + r"|\d+|.", re.S)


def bracket_end(code, i):
    """The position after the bracket that closes the one opened at i."""
    depth = 0
    for m in TOKEN.finditer(code, i):
        if m.group(0) in "[({":
            depth += 1
        elif m.group(0) in "])}":
            depth -= 1
            if depth == 0:
                return m.end()
    return len(code)


def colours(code, known):
    """The colour of each character of the code, and the names it gives a value to. `known`: the names with a value
    (the package's once it is loaded, those set by the cells before), with those of the system (SYSTEM_NAMES); the
    others are undefined; the variables of a Sum are local inside it, and the names of patterns (i_) in the line that
    has them."""
    local = []      # (start, end, names)
    for m in re.finditer(r"\bSum\[", code):
        end = bracket_end(code, m.end() - 1)
        local.append((m.start(), end, set(re.findall(r",\s*\{\s*(" + NAME + r")\s*,", code[m.end():end]))))
    start = 0
    for line in code.split("\n"):
        local.append((start, start + len(line), set(re.findall("(" + NAME + ")_", line))))
        start += len(line) + 1
    out = []
    for m in TOKEN.finditer(code):
        t = m.group(0)
        if t.startswith("(*"):
            kind = "comment"
        elif t.startswith('"'):
            kind = "string"
        elif re.match(NAME, t):
            if any(a <= m.start() < b and t in names for a, b, names in local):
                kind = "local"
            elif t in known or t in SYSTEM_NAMES:
                kind = "plain"
            else:
                kind = "symbol"
        else:
            kind = "plain"
        out.extend([CODE_COLORS[kind]] * len(t))
    assigned = set(re.findall(r"^\s*(" + NAME + r")\s*=(?!=)", code, re.M))
    return out, assigned


def expressions(code):
    """The number of expressions of a cell: a notebook evaluates each line that ends one on its own, with its own number
    (In[n]:= is that of the first one, Out[n]= that of the one that gives the output)."""
    flat = re.sub(r"\(\*.*?\*\)|\"(?:[^\"\\]|\\.)*\"", "0", code, flags=re.S)
    n = depth = 0
    for line in flat.split("\n"):
        depth += sum(line.count(b) for b in "[({") - sum(line.count(b) for b in "])}")
        if depth == 0 and line.strip():
            n += 1
    return n


def wrap(code, cols, indent=2):
    """The row and column of each character of the code (and of its end) in a cell `cols` characters wide: a line too
    long is wrapped as a notebook wraps it, where the nesting is least deep (after a comma, before a binary operator),
    the rest indented by `indent`."""
    pos, row = [], 0
    for line in code.split("\n"):
        start, ind = 0, 0
        while len(line) - start + ind > cols:
            best, depth, instr = None, 0, False
            for j in range(start, len(line)):
                ch = line[j]
                if ch == '"':
                    instr = not instr
                if instr:
                    continue
                if j > start + 1 and j - start + ind <= cols and line[j - 1] == " " and (
                        line[j - 2] == "," or (line[j] in "+-" and line[j + 1:j + 2] == " ") or line.startswith("**", j)):
                    if best is None or depth <= best[0]:
                        best = (depth, j)
                depth += (ch in "[({") - (ch in "])}")
            if best is None:
                break
            pos.extend((row, ind + k) for k in range(best[1] - start))
            row, start, ind = row + 1, best[1], indent
        pos.extend((row, ind + k) for k in range(len(line) - start + 1))      # with the newline (or the end)
        row += 1
    return pos


class CodeCell:
    """An input cell: the code, coloured, at (x, y) (the top of its first row), wrapped to the width of the cell."""

    def __init__(self, code, known, fnt, x, y, line_height, cols):
        self.code, self.font, self.x, self.y, self.lh = code, fnt, x, y, line_height
        self.colours, self.assigned = colours(code, known)
        self.pos = wrap(code, cols)
        self.advance = fnt.getlength("M")
        self.rows = self.pos[-1][0] + 1
        self.height = self.rows * line_height
        self.ascent = line_height * 0.74     # from the top of a row to the baseline

    def baseline(self, row):
        return self.y + row * self.lh + self.ascent

    def draw(self, draw, n, cursor=False):
        """Draw the first n characters; with the cursor after them."""
        run, last = "", None
        for i in range(n + 1):
            key = (self.pos[i][0], self.colours[i]) if i < n and self.code[i] != "\n" else None
            if key != last or (key and self.pos[i][1] != self.pos[i - 1][1] + 1):
                if run.strip():
                    row, col = self.pos[i - len(run)]
                    draw.text((self.x + col * self.advance, self.baseline(row)), run, font=self.font, fill=last[1], anchor="ls")
                run = ""
            if key:
                run += self.code[i]
            last = key
        if cursor:
            if n and self.code[n - 1] != "\n":
                row, col = self.pos[n - 1][0], self.pos[n - 1][1] + 1
            else:
                row, col = self.pos[n]
            cx, base = self.x + col * self.advance + px(1), self.baseline(row)
            draw.rectangle([cx, base - px(16), cx + ipx(2) - 1, base + px(5)], fill=ACCENT)


# ---------------------------------------------------------------------------------------------------------------------
# frames

def ease(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3 - 2 * u)


def ease_out(u):
    """Quick at first, then settling: the fade-in of the elements."""
    u = min(max(u, 0.0), 1.0)
    return 1 - (1 - u) ** 3


def ease_out_back(u, overshoot=1.6):
    """Like ease_out, a little beyond 1 before settling: the logo popping in."""
    u = min(max(u, 0.0), 1.0) - 1
    return 1 + (overshoot + 1) * u ** 3 + overshoot * u ** 2


def smootherstep(u):
    """0 to 1 with zero slope and zero curvature at both ends."""
    u = min(max(u, 0.0), 1.0)
    return u * u * u * (u * (6 * u - 15) + 10)


def bump(u):
    """0 at u = 0 and 1, up to 1 in the middle, with zero slope at both ends: a pulse."""
    return 16 * u * u * (1 - u) * (1 - u) if 0 < u < 1 else 0.0


def with_alpha(img, a):
    if a >= 1:
        return img
    out = img.copy()
    out.putalpha(img.getchannel("A").point(lambda v: int(v * a)))
    return out


class Item:
    """An image on the frame from time t0, fading in (and rising a few pixels), until t1 (or the end of its scene)."""

    def __init__(self, img, x, y, t0=0.0, fade=0.28, t1=None, rise=px(10)):
        self.img, self.x, self.y, self.t0, self.fade, self.t1, self.rise = img, x, y, t0, fade, t1, rise

    def draw(self, frame, t, dy=0):
        if t < self.t0 or (self.t1 is not None and t >= self.t1):
            return
        a = ease_out((t - self.t0) / self.fade) if self.fade else 1
        y = int(round(self.y - dy + (1 - a) * self.rise))
        if y > H or y + self.img.height < 0:
            return
        paste_clipped(frame, with_alpha(self.img, a), int(round(self.x)), y)


# ---------------------------------------------------------------------------------------------------------------------
# the logo, from Logo/FermiWick_logo.svg

def _colour(v):
    return tuple(int(v[i:i + 2], 16) for i in (1, 3, 5))


def _numbers(d):
    """The commands and numbers of the path data d (absolute commands only)."""
    return re.findall(r"[A-Za-z]|-?\d*\.?\d+", d)


class Bracket:
    """A contraction of the logo: its path (lines and quadratic corners) sampled densely, from the creator (the open
    circle) to the annihilator (the filled circle), with the arc length along it."""

    def __init__(self, d, colour, width):
        self.colour, self.width = colour, width
        toks, pts, i, cmd = _numbers(d), [], 0, None
        while i < len(toks):
            if toks[i].isalpha():
                cmd, i = toks[i], i + 1
            if cmd in "Zz":
                break
            if cmd == "M":
                pts, i = [(float(toks[i]), float(toks[i + 1]))], i + 2
            elif cmd in "LHV":
                x, y = pts[-1]
                if cmd == "L":
                    end, i = (float(toks[i]), float(toks[i + 1])), i + 2
                else:
                    end, i = ((float(toks[i]), y) if cmd == "H" else (x, float(toks[i]))), i + 1
                n = max(1, int(math.ceil(math.dist(pts[-1], end) / 0.25)))
                pts += [(x + (end[0] - x) * k / n, y + (end[1] - y) * k / n) for k in range(1, n + 1)]
            elif cmd == "Q":
                (x0, y0), (cx, cy), (x1, y1) = pts[-1], (float(toks[i]), float(toks[i + 1])), (float(toks[i + 2]), float(toks[i + 3]))
                i += 4
                pts += [((1 - u) ** 2 * x0 + 2 * u * (1 - u) * cx + u * u * x1, (1 - u) ** 2 * y0 + 2 * u * (1 - u) * cy + u * u * y1)
                        for u in (k / 96 for k in range(1, 97))]
            else:
                sys.exit("Logo/FermiWick_logo.svg: path command %s not read by this script" % cmd)
        self.points = pts
        self.orient()

    def orient(self, reverse=False):
        if reverse:
            self.points = self.points[::-1]
        self.lengths = [0.0]
        for a, b in zip(self.points, self.points[1:]):
            self.lengths.append(self.lengths[-1] + math.dist(a, b))
        self.length = self.lengths[-1]

    def point(self, s):
        """The point at the arc length s and the unit tangent there (the direction of travel)."""
        s = min(max(s, 0.0), self.length)
        lo, hi = 0, len(self.lengths) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.lengths[mid] <= s:
                lo = mid
            else:
                hi = mid
        (x0, y0), (x1, y1) = self.points[lo], self.points[hi]
        seg = self.lengths[hi] - self.lengths[lo]
        u = (s - self.lengths[lo]) / seg if seg else 0.0
        return x0 + u * (x1 - x0), y0 + u * (y1 - y0), (x1 - x0) / seg, (y1 - y0) / seg

    def stroke(self, s):
        """The outline of the stroke from the creator to the arc length s (butt ends), as a polygon."""
        left, right, hw = [], [], self.width / 2
        samples = [k for k in range(len(self.lengths)) if self.lengths[k] < s] + [None]
        for k in samples:
            if k is None:
                x, y, tx, ty = self.point(s)
            else:
                x, y = self.points[k]
                a, b = self.points[max(k - 1, 0)], self.points[min(k + 1, len(self.points) - 1)]
                n = math.dist(a, b)
                tx, ty = (b[0] - a[0]) / n, (b[1] - a[1]) / n
            left.append((x - ty * hw, y + tx * hw))
            right.append((x + ty * hw, y - tx * hw))
        return left + right[::-1]


class Head:
    """An arrow head of the logo (a triangle stroked with round joins), on its bracket at the arc length `rest`, in the
    frame of the bracket there (along it, across it), so that it can travel along the bracket, turning with it."""

    def __init__(self, d, colour, width, bracket):
        toks = _numbers(d)
        pts = [(float(toks[i]), float(toks[i + 1])) for i in (1, 4, 7)]       # M x y L x y L x y Z
        self.colour, self.radius, self.bracket = colour, width / 2, bracket
        # the tip: the vertex between the two equal sides
        sides = lambda k: abs(math.dist(pts[k], pts[k - 1]) - math.dist(pts[k], pts[(k + 1) % 3]))
        tip = min(range(3), key=sides)
        base = [pts[k] for k in range(3) if k != tip]
        mid = ((base[0][0] + base[1][0]) / 2, (base[0][1] + base[1][1]) / 2)
        centre = ((pts[tip][0] + mid[0]) / 2, (pts[tip][1] + mid[1]) / 2)
        self.rest = min((math.dist(centre, p), s) for p, s in zip(bracket.points, bracket.lengths))[1]
        x, y, tx, ty = bracket.point(self.rest)
        if (pts[tip][0] - mid[0]) * tx + (pts[tip][1] - mid[1]) * ty <= 0:
            sys.exit("Logo/FermiWick_logo.svg: an arrow head does not point from the creator to the annihilator")
        self.local = [((p[0] - x) * tx + (p[1] - y) * ty, -(p[0] - x) * ty + (p[1] - y) * tx) for p in pts]

    def polygon(self, s, g=1.0):
        """The triangle at the arc length s, scaled by g about its centre."""
        x, y, tx, ty = self.bracket.point(s)
        return [(x + g * (u * tx - v * ty), y + g * (u * ty + v * tx)) for u, v in self.local]


class Logo:
    """The geometry of Logo/FermiWick_logo.svg, in units of its 128 x 128 picture: the badge (a rounded square in a
    gradient, a soft white light, a faint white edge), the two contractions (Bracket), their arrow heads (Head) and the
    four operators (circles: filled for an annihilator, open for a creator)."""

    def __init__(self, path):
        root = ET.parse(path).getroot()
        tag = lambda e: e.tag.rsplit("}", 1)[-1]
        if root.get("viewBox") != "0 0 128 128":
            sys.exit("Logo/FermiWick_logo.svg: a viewBox of 128 x 128 expected")
        paths, heads, self.circles, rects = [], [], [], []
        for e in root.iter():
            t = tag(e)
            if t == "linearGradient":
                stops = [_colour(s.get("stop-color")) for s in e if tag(s) == "stop"]
                self.start, self.end = stops[0], stops[-1]
                self.axis = tuple(float(e.get(k)) for k in ("x1", "y1", "x2", "y2"))
            elif t == "radialGradient":
                stops = [s for s in e if tag(s) == "stop"]
                self.light = tuple(float(e.get(k)) for k in ("cx", "cy", "r")) + (float(stops[0].get("stop-opacity")),)
            elif t == "rect":
                rects.append(e)
            elif t == "path" and e.get("fill") not in (None, "none"):
                heads.append((e.get("d"), _colour(e.get("fill")), float(e.get("stroke-width", "0"))))
            elif t == "path":
                paths.append(Bracket(e.get("d"), _colour(e.get("stroke")), float(e.get("stroke-width"))))
            elif t == "circle":
                filled = e.get("fill") not in (None, "none")
                self.circles.append(dict(c=(float(e.get("cx")), float(e.get("cy"))), r=float(e.get("r")), filled=filled,
                                         colour=_colour(e.get("fill") if filled else e.get("stroke")),
                                         width=0.0 if filled else float(e.get("stroke-width"))))
        self.radius = float(rects[0].get("rx"))
        edge = next(r for r in rects if r.get("stroke"))
        self.edge = (float(edge.get("x")), float(edge.get("rx")), float(edge.get("stroke-width")), float(edge.get("stroke-opacity")))
        # each bracket from its creator to its annihilator: the circles of its colour at its two ends
        self.brackets = paths
        for b in paths:
            ends = [min((c for c in self.circles if c["colour"] == b.colour), key=lambda c: math.dist(c["c"], p))
                    for p in (b.points[0], b.points[-1])]
            if ends[0]["filled"] == ends[1]["filled"]:
                sys.exit("Logo/FermiWick_logo.svg: a contraction does not join a creator to an annihilator")
            b.orient(reverse=ends[0]["filled"])
            b.creator, b.annihilator = (ends[1], ends[0]) if ends[0]["filled"] else (ends[0], ends[1])
        self.heads = [Head(d, colour, width, next(b for b in paths if b.colour == colour)) for d, colour, width in heads]


LOGO = Logo(os.path.join(ROOT, "Logo", "FermiWick_logo.svg"))
LOGO_SUPERSAMPLE = 4     # the logo is drawn this many times larger and reduced, for smooth edges


@functools.lru_cache(maxsize=64)
def badge_image(big):
    """The badge of the logo at `big` px, without the drawing: the rounded square in the gradient of the SVG (its
    axis from the top left to the bottom right), under the soft white light of its radial gradient (falling linearly to
    none), and the faint white edge. Both gradients are smooth: they are computed on a grid of 256 x 256 (two points per
    unit) and enlarged. The badge does not move: one image per size, the same in every frame."""
    f = big / 128.0
    n = 256
    units = [(i + 0.5) * 128.0 / n for i in range(n)]
    x1, y1, x2, y2 = LOGO.axis
    ax, ay = x2 - x1, y2 - y1
    along = bytes(min(255, max(0, int(255 * ((x - x1) * ax + (y - y1) * ay) / (ax * ax + ay * ay) + 0.5)))
                  for y in units for x in units)
    cx, cy, r, opacity = LOGO.light
    light = bytes(int(255 * opacity * max(0.0, 1 - math.hypot(x - cx, y - cy) / r) + 0.5) for y in units for x in units)
    along = Image.frombytes("L", (n, n), along).resize((big, big), Image.BILINEAR)
    light = Image.frombytes("L", (n, n), light).resize((big, big), Image.BILINEAR)
    white = Image.new("L", (big, big), 255)
    badge = Image.composite(Image.new("RGBA", (big, big), LOGO.end + (255,)),
                            Image.new("RGBA", (big, big), LOGO.start + (255,)), along)
    badge.alpha_composite(Image.merge("RGBA", (white, white, white, light)))
    # the edge: a stroke centred on a rounded square inset by x (rx), here from the outline of the badge inwards
    x, rx, width, opacity = LOGO.edge
    edge = Image.new("L", (big, big), 0)
    ImageDraw.Draw(edge).rounded_rectangle([(x - width / 2) * f, (x - width / 2) * f, big - 1 - (x - width / 2) * f,
                                            big - 1 - (x - width / 2) * f], radius=int(round((rx + width / 2) * f)),
                                           outline=round(opacity * 255), width=int(round(width * f)))
    badge.alpha_composite(Image.merge("RGBA", (white, white, white, edge)))
    shape = Image.new("L", (big, big), 0)
    ImageDraw.Draw(shape).rounded_rectangle([0, 0, big - 1, big - 1], radius=int(round(LOGO.radius * f)), fill=255)
    badge.putalpha(shape)
    return badge


class LogoState:
    """What the logo shows: how much of each bracket is drawn (arc length from its creator), where each arrow head is
    (arc length) and how large (1 in place), how large each circle is (1 in place). At rest, the logo of the SVG."""

    def __init__(self, fronts=None, heads=None, circles=None):
        self.fronts = fronts or [b.length for b in LOGO.brackets]
        self.heads = heads or [(h.rest, 1.0) for h in LOGO.heads]
        self.circles = circles or [1.0] * len(LOGO.circles)

    def key(self):
        return tuple(round(v, 4) for v in self.fronts + [v for h in self.heads for v in h] + self.circles)


REST = LogoState().key()


@functools.lru_cache(maxsize=8)
def rest_logo(size):
    return draw_logo(size, LogoState())


def logo_image(size, state):
    return rest_logo(size) if state.key() == REST else draw_logo(size, state)


def draw_logo(size, state):
    """The logo at `size` px in the state `state`, drawn LOGO_SUPERSAMPLE times larger and reduced. Its strokes are
    painted through one mask for each colour (the brackets, heads and circles of one colour never touch those of the
    other), so that overlapping shapes of one colour make one shape."""
    big = size * LOGO_SUPERSAMPLE
    f = big / 128.0
    img = badge_image(big).copy()
    masks = {}

    def T(x, y):
        # PIL puts the centre of a pixel at its integer coordinates; SVG at half-integers
        return x * f - 0.5, y * f - 0.5

    def paint(colour):
        if colour not in masks:
            masks[colour] = Image.new("L", (big, big), 0)
        return ImageDraw.Draw(masks[colour])

    def disc(d, x, y, r):
        (X, Y), R = T(x, y), r * f
        d.ellipse([X - R, Y - R, X + R, Y + R], fill=255)

    for b, front in zip(LOGO.brackets, state.fronts):
        if front > 0.05:
            paint(b.colour).polygon([T(x, y) for x, y in b.stroke(front)], fill=255)
    for h, (s, g) in zip(LOGO.heads, state.heads):
        if g > 0.02:
            d, pts, r = paint(h.colour), h.polygon(s, g), h.radius * g
            d.polygon([T(x, y) for x, y in pts], fill=255)
            for k in range(3):
                (xa, ya), (xb, yb) = pts[k], pts[(k + 1) % 3]
                n = math.dist(pts[k], pts[(k + 1) % 3])
                nx, ny = -(yb - ya) / n * r, (xb - xa) / n * r
                d.polygon([T(xa + nx, ya + ny), T(xb + nx, yb + ny), T(xb - nx, yb - ny), T(xa - nx, ya - ny)], fill=255)
                disc(d, xa, ya, r)
    for c, k in zip(LOGO.circles, state.circles):
        if k > 0.02:
            d = paint(c["colour"])
            if c["filled"]:
                disc(d, c["c"][0], c["c"][1], c["r"] * k)
            else:
                # a ring: PIL draws the outline inside the box, as wide as the stroke
                (X, Y), R = T(*c["c"]), (c["r"] + c["width"] / 2) * k * f
                d.ellipse([X - R, Y - R, X + R, Y + R], outline=255, width=max(1, int(round(c["width"] * k * f))))
    for colour, mask in masks.items():
        img.paste(colour + (255,), (0, 0, big, big), mask)
    return img.resize((size, size), Image.LANCZOS)


def emerging(s, length, reach=9.0):
    """The size of an arrow head at the arc length s of a bracket of `length`: it grows out of the creator and shrinks
    into the annihilator within `reach` units of them, 1 everywhere else."""
    return smootherstep(s / reach) * smootherstep((length - s) / reach)


class Turn:
    """The heads of the logo travelling: by `amount` (in lengths of their brackets: 1, so that they land on their rest
    positions), starting at the time `start` (s). The motion is that of a heavy wheel pulled by a soft spring: the pull
    runs with a velocity that rises from rest (smootherstep, so with no jerk), glides and comes back to rest; the heads
    follow it through a damped spring (natural frequency `omega`, damping ratio `zeta`), so that they start a little
    late and softly, glide, and land just beyond their rest positions before settling on them. The motion is integrated
    once into a table; the last fraction of the settling is faded out, so that the heads end exactly at rest."""

    def __init__(self, start, amount=1.0, ramp_up=0.5, glide=0.45, ramp_down=0.5, omega=10.0, zeta=0.5, settle=0.55):
        self.start, self.amount = start, amount
        self.end = ramp_up + glide + ramp_down          # the end of the pull, from the start
        self.length = self.end + settle                 # at rest from then on
        self.dt = 1 / 2400
        top = amount / (ramp_up / 2 + glide + ramp_down / 2)

        def pull(t):
            if t < ramp_up:
                return top * smootherstep(t / ramp_up)
            if t < ramp_up + glide:
                return top
            return top * (1 - smootherstep((t - ramp_up - glide) / ramp_down))

        x = th = w = 0.0
        self.table = [0.0]
        for k in range(int(self.length / self.dt) + 2):
            t = k * self.dt
            # semi-implicit Euler: the pull x, the heads th (velocity w)
            x += pull(t) * self.dt
            w += (omega ** 2 * (x - th) - 2 * zeta * omega * w) * self.dt
            th += w * self.dt
            fade = smootherstep((t - self.end) / settle)       # the residual of the settling, faded out
            self.table.append(amount + (th - amount) * (1 - fade))

    def __call__(self, t):
        """How far the heads have travelled at the time t (in lengths of their brackets)."""
        u = (t - self.start) / self.dt
        if u <= 0:
            return self.table[0]
        k = int(u)
        if k >= len(self.table) - 1:
            return self.table[-1]
        return self.table[k] + (u - k) * (self.table[k + 1] - self.table[k])

    def when(self, value):
        """The first time the heads have travelled `value`."""
        return self.start + next(k for k, v in enumerate(self.table) if v >= value) * self.dt


PULSE, PULSE_TIME = 0.16, 0.34      # how much a circle swells when a head goes in or out of it, and for how long


class Flow:
    """The title card: the logo at rest, then from `start` its arrow heads travel once along their contractions (Turn),
    from their places into the annihilators, out of the creators again and back to their places; each circle swells
    a little as a head goes into it or out of it."""

    def __init__(self, start):
        self.turn = Turn(start)
        self.pulses = [[] for _ in LOGO.circles]
        for h in LOGO.heads:
            b = h.bracket
            when = self.turn.when(1 - h.rest / b.length)
            for c in (b.creator, b.annihilator):
                self.pulses[LOGO.circles.index(c)].append(when)

    def state(self, t):
        u, heads = self.turn(t), []
        for h in LOGO.heads:
            length = h.bracket.length
            s = (h.rest + u * length) % length
            heads.append((s, emerging(s, length)))
        circles = [1 + PULSE * sum(bump((t - w) / PULSE_TIME) for w in ws) for ws in self.pulses]
        return LogoState(heads=heads, circles=circles)


class DrawIn:
    """The closing card: the logo draws itself. From `start` the four circles appear from left to right (the operators,
    each popping a little beyond its size and back), then each contraction is drawn from its creator to its annihilator
    (the pen with a soft start and stop), its arrow head riding behind the pen on a stiff spring and settling on its
    place as the pen runs on; the annihilator swells a little when the pen reaches it."""

    def __init__(self, start, stagger=0.08, pop=0.36, pen=0.62, draw=0.8, omega=17.0, zeta=0.55):
        self.start, self.pop, self.pen, self.draw = start, pop, pen, draw
        order = sorted(range(len(LOGO.circles)), key=lambda i: LOGO.circles[i]["c"][0])
        self.appear = {i: start + stagger * k for k, i in enumerate(order)}
        self.arrive = start + pen + draw
        self.dt, self.tables = 1 / 2400, []
        settle = 0.5
        for h in LOGO.heads:
            x = w = 0.0
            table = [0.0]
            for k in range(int((pen + draw + settle) / self.dt) + 2):
                t = k * self.dt
                target = min(self.front(start + t) * h.bracket.length, h.rest)
                w += (omega ** 2 * (target - x) - 2 * zeta * omega * w) * self.dt
                x += w * self.dt
                fade = smootherstep((t - pen - draw) / settle)      # the residual of the settling, faded out
                table.append(max(0.0, h.rest + (x - h.rest) * (1 - fade)))
            self.tables.append(table)

    def front(self, t):
        """The fraction of each bracket drawn at the time t."""
        return smootherstep((t - self.start - self.pen) / self.draw)

    def state(self, t):
        p = self.front(t)
        heads = []
        for h, table in zip(LOGO.heads, self.tables):
            k = min(max(int((t - self.start) / self.dt), 0), len(table) - 1)
            heads.append((table[k], emerging(table[k], h.bracket.length)))
        annihilators = [LOGO.circles.index(b.annihilator) for b in LOGO.brackets]
        circles = [ease_out_back((t - self.appear[i]) / self.pop, 1.2) if t < self.appear[i] + self.pop else 1.0
                   for i in range(len(LOGO.circles))]
        for i in annihilators:
            circles[i] += PULSE * bump((t - self.arrive + 0.04) / PULSE_TIME)
        return LogoState(fronts=[p * b.length for b in LOGO.brackets], heads=heads,
                         circles=[max(0.0, c) for c in circles])


class LogoItem:
    """The logo, alive: it pops in (scaled up from 80%, a little beyond its size and back) at t0, then animates
    (`anim`, a Flow or a DrawIn: its state at each time)."""

    def __init__(self, size, x, y, t0=0.0, fade=0.45, anim=None):
        self.size, self.x, self.y, self.t0, self.fade, self.anim = size, x, y, t0, fade, anim

    def draw(self, frame, t, dy=0):
        if t < self.t0:
            return
        u = (t - self.t0) / self.fade if self.fade else 1.0
        scale = 0.8 + 0.2 * ease_out_back(u)
        size = max(ipx(8), int(round(self.size * scale)))
        img = with_alpha(logo_image(size, self.anim.state(t) if self.anim else LogoState()), ease_out(1.6 * u))
        paste_clipped(frame, img, int(round(self.x + (self.size - size) / 2)), int(round(self.y - dy + (self.size - size) / 2)))


def paste_clipped(frame, img, x, y):
    top, bottom = max(0, -y), min(img.height, H - y)
    if bottom > top:
        frame.alpha_composite(img.crop((0, top, img.width, bottom)), (x, y + top))


def blank():
    frame = Image.new("RGBA", (W, H), WHITE + (255,))
    ImageDraw.Draw(frame).rectangle([0, H - ipx(6), W, H], fill=ACCENT)   # the band at the foot
    return frame


def fade_to_white(frame, a):
    """a = 1: the frame; a = 0: white (the accent band stays)."""
    if a >= 1:
        return frame
    return Image.blend(blank(), frame, a)


class Scene:
    def __init__(self, duration, items=(), draw=None, fade_out=0.22, full=0.0):
        self.duration, self.items, self.custom, self.fade_out = duration, list(items), draw, fade_out
        self.full = full        # from this time on, everything of the scene is in view

    def frame(self, t):
        frame = blank()
        if self.custom:
            self.custom(frame, t)
        for item in self.items:
            item.draw(frame, t)
        if self.fade_out and t > self.duration - self.fade_out:
            a = (self.duration - t) / self.fade_out
            shift = int(round(px(36) * (1 - a) ** 2))
            if shift:
                moved = blank()
                moved.paste(frame.crop((shift, 0, W, H - ipx(6))), (0, 0))
                frame = moved
            frame = fade_to_white(frame, a)
        return frame


def centred(img, y):
    return (W - img.width) / 2, y


def title_scene(el, duration, parts, anim, fade_out=0.3, shown=False):
    """The logo and, under it, the elements `parts` (the name, the tagline or the address of the site), stacked in the
    middle; the logo animated by `anim` (Flow, DrawIn). With `shown`, the card is in full view from its first frame (the
    first frame of the video, which a player shows before it plays) instead of building up."""
    size = round(2.6 / 2.54 * DPI)
    slot = size + 2 * BORDER * PT                 # the logo as tall as an element 2.6 cm high
    gaps = {"name": px(26), "tagline": px(12), "url": px(30)}
    height = slot + sum(el[p].height + gaps[p] for p in parts)
    y = (H - height) / 2 - px(8)
    times = {"name": 0.25, "tagline": 0.45, "url": 0.7}
    items = [LogoItem(size, (W - size) / 2, y + (slot - size) / 2, 0.0, fade=0 if shown else 0.45, anim=anim)]
    y += slot
    for p in parts:
        y += gaps[p]
        x, _ = centred(el[p], y)
        items.append(Item(el[p], x, y, 0.0 if shown else times[p], fade=0 if shown else 0.35))
        y += el[p].height
    return Scene(duration, items, fade_out=fade_out)


def caption_lines(el, base):
    """The baselines of the headings and of their second lines, the same on every card (the top of the tallest heading
    44 px from the top, as in the cards of FunctionalDiagrams), and the top of the space left under them."""
    heads = [k for k in el if re.match(r"h\d+$", k)]
    subs = [k for k in el if re.match(r"s\d+$", k)]
    head = px(44) + max(base[k] for k in heads)
    sub = head + max(el[k].height - base[k] for k in heads) + px(2) + max(base[k] for k in subs)
    return head, sub, sub + max(el[k].height - base[k] for k in subs) + px(30)


class Card:
    """One card of the session: its caption (top left: the heading and its second line), then its cells typed one after
    the other and evaluated, each followed by its output, as in a notebook: In[n]:= on the left of an input cell once it
    is evaluated, Out[n]= on the left of the output, and on the right a bracket for each cell and one around the group
    of an input with its output. The cells and their outputs are left-aligned with the caption, a little above the
    middle of the height left under it."""
    X_CODE, LH, CELL_GAP, OUT_GAP = px(122), px(26), px(16), px(12)
    X_BRACKET, X_GROUP = W - ipx(36), W - ipx(27)     # the cell brackets, and those of the groups

    def __init__(self, el, base, k, cells, session, hold, start=0.45):
        """cells: [(code, name of the output element or None)]; session: the number of the next input line and the
        names with a value, carried from card to card."""
        self.codefont = font("SourceCodePro-Medium.otf", CODE_SIZE)
        self.labelfont = font("SourceCodePro-Regular.otf", px(14))
        cols = int((self.X_BRACKET - px(16) - self.X_CODE) / self.codefont.getlength("M"))
        border = BORDER * PT
        head, sub, top = caption_lines(el, base)
        h, s = "h%d" % k, "s%d" % k
        self.items = [Item(el[h], self.X_CODE - border, head - base[h], 0.0, fade=0.22),
                      Item(el[s], self.X_CODE - border, sub - base[s], 0.07, fade=0.22)]
        for e in (h, s):
            if self.X_CODE + el[e].width - 2 * border > self.X_BRACKET - px(16):
                print("warning: the caption %s of card %d is too wide" % (e, k))
        # the cells, laid out from 0 and then moved down to their place
        self.cells, t, y, self.full = [], start, 0.0, start
        for code, out in cells:
            cell = CodeCell(code, session["known"], self.codefont, self.X_CODE, y, self.LH, cols)
            n = expressions(code)
            typing = min(0.95, max(0.35, len(code) / 220))
            c = dict(cell=cell, n=session["line"], t0=t, cps=len(code) / typing, t_eval=t + typing + 0.12, top=y - px(5))
            session["line"] += n
            session["known"] |= cell.assigned
            if re.search(r"Needs\[\"FermiWick`\"\]", code):
                # the names of the package have a value from the next cell on
                session["known"] |= package_names()
            y += cell.height
            c["input_bottom"] = y + px(2)
            t = c["t_eval"] + 0.18
            c["out"] = None
            if out:
                y += self.OUT_GAP
                img = el[out]
                if self.X_CODE + img.width - 2 * border > self.X_BRACKET - px(16):
                    print("warning: the output %s runs into the cell brackets" % out)
                # Out[n]= with the number of the last expression, on the baseline of the output
                c["out"] = (Item(img, self.X_CODE - border, y, t, fade=0.3), session["line"] - 1, base[out])
                self.full = t + 0.3
                y += img.height - 2 * border
            t += 0.22
            c["bottom"] = y + px(4)
            y += self.CELL_GAP
            self.cells.append(c)
        height = y - self.CELL_GAP
        if height > H - px(40) - top:
            print("warning: the cells of card %d are %.0f px too tall" % (k, (height - (H - px(40) - top)) / S))
        dy = top + max(0.0, (H - px(40) - top - height) * 0.32)
        for c in self.cells:
            c["cell"].y += dy
            for key in ("top", "input_bottom", "bottom"):
                c[key] += dy
            if c["out"]:
                c["out"][0].y += dy - border
        self.height = height
        self.duration = self.full + hold

    def draw(self, frame, t):
        d = ImageDraw.Draw(frame)
        for item in self.items:
            item.draw(frame, t)
        for c in self.cells:
            if t < c["t0"]:
                break
            cell = c["cell"]
            n = min(len(cell.code), int((t - c["t0"]) * c["cps"]))
            typing = t < c["t_eval"]
            blink = n < len(cell.code) or int(t * 2.2) % 2 == 0
            cell.draw(d, n, cursor=typing and blink)
            if t >= c["t_eval"]:
                self.label(d, "In[%d]:=" % c["n"], cell.baseline(0), LABEL)
            self.bracket(d, self.X_BRACKET, c["top"], c["input_bottom"])
            if c["out"] and t >= c["out"][0].t0:
                item, number, base = c["out"]
                item.draw(frame, t)
                a = ease((t - item.t0) / 0.3)
                colour = tuple(int(255 - (255 - v) * a) for v in LABEL)
                rise = (1 - ease_out((t - item.t0) / item.fade)) * item.rise
                self.label(d, "Out[%d]=" % number, item.y + base + rise, colour)
                self.bracket(d, self.X_BRACKET, item.y + BORDER * PT - px(5), c["bottom"],
                             colour=tuple(int(255 - (255 - v) * a) for v in BRACKET))
                self.bracket(d, self.X_GROUP, c["top"], c["bottom"])

    def label(self, d, text, baseline, colour):
        """In[n]:= or Out[n]=, right-aligned before the cells, on the baseline given."""
        d.text((self.X_CODE - px(12), baseline), text, font=self.labelfont, fill=colour, anchor="rs")

    @staticmethod
    def bracket(d, x, top, bottom, colour=BRACKET):
        d.line([x, top, x, bottom], fill=colour, width=ipx(1))
        d.line([x - ipx(5), top, x, top], fill=colour, width=ipx(1))
        d.line([x - ipx(5), bottom, x, bottom], fill=colour, width=ipx(1))


# ---------------------------------------------------------------------------------------------------------------------

def timeline(el, base, session_cards, preview=False):
    """The scenes: the title card, the cards of the session (CARDS), the closing card. The preview has no title card,
    shorter holds, and its closing card fades out (it loops)."""
    session = {"line": 1, "known": set()}
    cards = []
    for k, (_, _, cells, hold) in enumerate(session_cards, 1):
        cells = [(code, output_name(k, j) if out else None) for j, (code, out) in enumerate(cells)]
        card = Card(el, base, k, cells, session, hold - (0.6 if preview else 0))
        cards.append(Scene(card.duration, draw=card.draw, full=card.full))
    # the closing card is the title card, so that the last frame of the video is its first and the loop has no seam
    end = title_scene(el, 2.8 if preview else 3.4, ["name", "tagline", "url"], DrawIn(0.2), fade_out=0.22 if preview else 0)
    if preview:
        return cards + [end]
    # the title card in full view from the first frame, which a player shows before the video plays
    title = title_scene(el, 2.6, ["name", "tagline", "url"], Flow(0.35), shown=True)
    return [title] + cards + [end]


def preview_start(first):
    """The preview loops from this time (s) of its first card, when all of it is in view: the first frame of the
    preview, the one seen by visitors who do not play animated images and in link previews, is that card (not the white
    frame before it builds up); its build-up, after the closing card, closes the loop on that frame."""
    return first.full + 0.2


_JOB = None    # (scenes, size) of the frames being drawn, read by the worker processes (inherited by fork)


def _frame_bytes(job):
    i, t = job
    scenes, size = _JOB
    f = scenes[i].frame(t).convert("RGB")
    return (f if size is None else f.resize(size, Image.LANCZOS)).tobytes()


def frames(scenes, fps, size=None, start=0.0):
    """The frames of the scenes one after the other, as raw RGB bytes; with `start` (s, in the first scene), a loop shown
    from there: the frames from that time to the end, then those before it. Each frame depends only on its time: they
    are drawn by several processes (fork) at once, in order, where fork is available."""
    global _JOB
    times = [(i, k / fps) for i, scene in enumerate(scenes) for k in range(int(round(scene.duration * fps)))]
    k0 = int(round(start * fps))
    times = times[k0:] + times[:k0]
    _JOB = (scenes, size)
    if "fork" not in multiprocessing.get_all_start_methods():
        yield from map(_frame_bytes, times)
        return
    with multiprocessing.get_context("fork").Pool(os.cpu_count() or 1) as pool:
        yield from pool.imap(_frame_bytes, times, chunksize=4)


def encode(scenes, fps, args, size=None, start=0.0):
    w, h = size or (W, H)
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "%dx%d" % (w, h),
           "-r", str(fps), "-i", "-"] + args
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for f in frames(scenes, fps, size, start):
        proc.stdin.write(f)
    proc.stdin.close()
    if proc.wait():
        sys.exit("ffmpeg failed")


def write_preview(scenes, path):
    """The looping preview, from preview_start, as an animated WebP at PREVIEW_FPS (written by PIL): a run of identical
    frames (a card in full view) is one frame, held as long; each frame is stored lossy or lossless, whichever is smaller
    (allow_mixed), so that the fades stay small and the lines of the cards sharp."""
    size, kept, counts = (PREVIEW_W, PREVIEW_H), [], []
    for f in frames(scenes, PREVIEW_FPS, size, preview_start(scenes[0])):
        if kept and f == kept[-1]:
            counts[-1] += 1
        else:
            kept.append(f)
            counts.append(1)
    # the delays in milliseconds, rounded so that they add up to the time of the frames (16 or 17 ms at 60 fps)
    durations, k = [], 0
    for c in counts:
        durations.append(round((k + c) * 1000 / PREVIEW_FPS) - round(k * 1000 / PREVIEW_FPS))
        k += c
    images = [Image.frombytes("RGB", size, f) for f in kept]
    images[0].save(path, save_all=True, append_images=images[1:], duration=durations, loop=0, minimize_size=True,
                   allow_mixed=True, quality=PREVIEW_QUALITY, method=6)


def write_poster(video):
    # the poster: the title card, the first frame of the video, so that nothing jumps when it starts (256 colours: the
    # gradient of the badge under its light shows bands with fewer)
    video[0].frame(0).convert("RGB").quantize(colors=256, method=Image.MEDIANCUT).save(
        os.path.join(MEDIA, "demo_poster.png"), optimize=True)


def main():
    os.makedirs(BUILD, exist_ok=True)
    if "--check" in sys.argv:
        check_session()
    os.makedirs(MEDIA, exist_ok=True)
    cards = session_cells()
    el, base = render_elements(cards)
    video = timeline(el, base, cards)
    total = sum(s.duration for s in video)
    print("video: %.1f s in %d scenes" % (total, len(video)))
    if "--poster" in sys.argv:
        write_poster(video)
        return
    preview = timeline(el, base, cards, preview=True)
    if "--frames" in sys.argv:
        out = os.path.join(BUILD, "frames")
        os.makedirs(out, exist_ok=True)
        for i, s in enumerate(video):
            for u in (0.55, 0.9):
                s.frame(u * s.duration).convert("RGB").save(os.path.join(out, "scene%d_%02d.png" % (i, round(u * 100))))
        for i, s in enumerate(preview):
            s.frame(0.9 * s.duration).convert("RGB").resize((PREVIEW_W, PREVIEW_H), Image.LANCZOS).save(os.path.join(out, "preview%d.png" % i))
        preview[0].frame(preview_start(preview[0])).convert("RGB").resize((PREVIEW_W, PREVIEW_H), Image.LANCZOS).save(
            os.path.join(out, "preview_first.png"))
        video[0].frame(0).convert("RGB").save(os.path.join(out, "scene0_first.png"))
        print("frames in", out)
        if "--only-frames" in sys.argv:
            return
    mp4 = os.path.join(MEDIA, "demo.mp4")
    encode(video, FPS, ["-c:v", "libx264", "-preset", "veryslow", "-crf", "20", "-tune", "animation", "-pix_fmt", "yuv420p",
                        "-movflags", "+faststart", mp4])
    write_poster(video)
    first = preview[0]
    if not preview_start(first) < first.duration - first.fade_out:
        print("warning: the first frame of the preview is not its first card in full view (preview_start)")
    webp = os.path.join(MEDIA, "demo_preview.webp")
    write_preview(preview, webp)
    for path, length in ((mp4, total), (webp, sum(s.duration for s in preview)), (os.path.join(MEDIA, "demo_poster.png"), 0)):
        print("%s: %.0f kB%s" % (os.path.relpath(path, ROOT), os.path.getsize(path) / 1024, ", %.1f s" % length if length else ""))


if __name__ == "__main__":
    main()
