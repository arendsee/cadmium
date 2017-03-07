TARGET=report

${TARGET}.pdf: ${TARGET}.Rnw R/*R
	Rscript -e "library(knitr); knit('${TARGET}.Rnw')"
	latexmk -pdf ${TARGET}.tex

.PHONY: tags
tags:
	ctags .

.PHONY: clean
clean:
	latexmk -f -c ${TARGET}.tex
	rm -f ${TARGET}.tex

.PHONY: deepclean
deepclean:
	latexmk -f -c ${TARGET}.tex 2> /dev/null
	rm -rf cache figure ${TARGET}.tex ${TARGET}.pdf *log  log.Rout
