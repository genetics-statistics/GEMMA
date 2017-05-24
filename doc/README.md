### Instructions for building PDF of GEMMA manual

The following commands will generate a PDF of the GEMMA manual from
the Latex:

```bash
pdflatex manual
bibtex manual
pdflatex manual
pdflatex manual
```

To run these commands, you will need a TeX distribution such as
[TexLive](https://www.tug.org/texlive) that includes commands
`pdflatex` and `bibtex`.

