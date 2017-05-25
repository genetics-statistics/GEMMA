# A Makefile for generating the PDF from Latex. Simply type 'make' in
# this directory to compile the Latex into a single PDF. Type 'make
# clean' to discard the output files except for the generated PDF.

# RULES
# -----
all: manual.pdf

# Create the PDF file.
manual.pdf : manual.tex manual.bib
	pdflatex manual
	bibtex manual
	pdflatex manual
	pdflatex manual

clean:
	rm -f manual.aux manual.bbl manual.blg
	rm -f manual.log manual.out manual.toc

