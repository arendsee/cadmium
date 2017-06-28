# NOTE: This is not a very reliable Makefile. It is convenient for me, but
# needs a major overhaul before general application.

TARGET=report.pdf
INPUT=input

SPECIES=${INPUT}/species

PROLOGUE=src/prologue
REPORT=src/report
EPILOGUE=src/epilogue

${TARGET}: species

.PHONY: load
load:
	cd ${PROLOGUE} && ./1_prepare-inputs.sh
	cd ${PROLOGUE} && ./2_extract-fasta.sh
	cd ${PROLOGUE} && ./3_gather-genome-summary-data.sh
	cd ${PROLOGUE} && ./4_get-search-intervals.sh

.PHONY: run
run:
	rm -rf cache
	cd ${REPORT} && ${MAKE} deepclean
	cd ${REPORT} && ${MAKE}

.PHONY: test
test:
	echo "Not yet implemented"

.PHONY: archive
archive:
	cd ${EPILOGUE} && ./archive.sh

.PHONY: clean
clean:
	rm -f *~
	rm -f *log log ${TARGET}
	rm -rf cache
	rm -rf ${INPUT}
	rm -f gmon.out
	rm -f preconfig.sh
	rm -f runconfig.R
	rm -f ${PROLOGUE}/config
	rm -f ${REPORT}/config
	rm -rf bin
	rm -rf src/smof src/synder
