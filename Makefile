# Makefile for article compilation
# Usage: make        -> compile English PDF
#        make fr     -> compile French PDF (main_fr.tex)
#        make clean  -> remove build artifacts (both languages)
#        make watch  -> continuous compilation (requires latexmk)

MAIN = main
MAIN_FR = main_fr

.PHONY: all fr clean watch view

all: $(MAIN).pdf

fr: $(MAIN_FR).pdf

$(MAIN).pdf: $(MAIN).tex references.bib
	pdflatex -interaction=nonstopmode $(MAIN).tex
	bibtex $(MAIN) || true
	pdflatex -interaction=nonstopmode $(MAIN).tex
	pdflatex -interaction=nonstopmode $(MAIN).tex

$(MAIN_FR).pdf: $(MAIN_FR).tex references.bib
	pdflatex -interaction=nonstopmode $(MAIN_FR).tex
	bibtex $(MAIN_FR) || true
	pdflatex -interaction=nonstopmode $(MAIN_FR).tex
	pdflatex -interaction=nonstopmode $(MAIN_FR).tex

clean:
	rm -f $(MAIN).pdf $(MAIN).aux $(MAIN).log $(MAIN).out \
		$(MAIN).toc $(MAIN).bbl $(MAIN).blg $(MAIN).nav $(MAIN).snm \
		$(MAIN).fls $(MAIN).fdb_latexmk $(MAIN).synctex.gz \
		$(MAIN_FR).pdf $(MAIN_FR).aux $(MAIN_FR).log $(MAIN_FR).out \
		$(MAIN_FR).toc $(MAIN_FR).bbl $(MAIN_FR).blg $(MAIN_FR).nav $(MAIN_FR).snm \
		$(MAIN_FR).fls $(MAIN_FR).fdb_latexmk $(MAIN_FR).synctex.gz

watch:
	latexmk -pdf -pvc -interaction=nonstopmode $(MAIN).tex

view: $(MAIN).pdf
	xdg-open $(MAIN).pdf 2>/dev/null || open $(MAIN).pdf
