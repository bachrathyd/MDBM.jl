# MDBM.jl LaTeX Documentation (Overleaf Compatible)

This folder contains the source for the MDBM.jl technical documentation.

## Files
- `main.tex`: The main LaTeX source file. It uses standard packages compatible with TeX Live and Overleaf.

## Usage on Overleaf
1. Create a new project on [Overleaf](https://www.overleaf.com).
2. Upload `main.tex` (and any image assets if you add them).
3. Set the "Main document" to `main.tex`.
4. Compile using `pdfLaTeX`.

## Local Compilation
To generate a PDF locally:

```bash
pdflatex main.tex
pdflatex main.tex # Run twice for Table of Contents
```

## Document Structure
The documentation covers:
- Basic and advanced API usage.
- Handling of periodic boundaries.
- Sub-cube interpolation techniques.
- Vector-valued parameter spaces.
- Practical engineering examples.