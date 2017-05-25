### Instructions for building PDF of GEMMA manual

The following commands will generate a PDF of the GEMMA manual from
the Latex:

```bash
cd doc
pdflatex manual
bibtex manual
pdflatex manual
pdflatex manual
```

Alternatively, you can use the Makefile:

```bash
cd doc
make
```

To run these commands, you will need a TeX distribution such as
[TexLive](https://www.tug.org/texlive) that includes commands
`pdflatex` and `bibtex`.
