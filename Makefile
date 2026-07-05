# Makefile for article compilation
# Usage: make        -> compile PDF
#        make clean  -> remove build artifacts
#        make watch  -> continuous compilation (requires latexmk)

MAIN = main

.PHONY: all clean watch view

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex references.bib
	pdflatex -interaction=nonstopmode $(MAIN).tex
	bibtex $(MAIN) || true
	pdflatex -interaction=nonstopmode $(MAIN).tex
	pdflatex -interaction=nonstopmode $(MAIN).tex

clean:
	rm -f $(MAIN).pdf $(MAIN).aux $(MAIN).log $(MAIN).out \
		$(MAIN).toc $(MAIN).bbl $(MAIN).blg $(MAIN).nav $(MAIN).snm \
		$(MAIN).fls $(MAIN).fdb_latexmk $(MAIN).synctex.gz

watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(MAIN).tex

view: $(MAIN).pdf
	xdg-open $(MAIN).pdf 2>/dev/null || open $(MAIN).pdf
