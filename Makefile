TARGET=fagin.loc

all:
	loc -kx locout ${TARGET}

.PHONY: run
run:
	./manifold-nexus.py main

.PHONY: clean
clean:
	rm -rf locout archive input manifold-nexus.py
	rm -f log
