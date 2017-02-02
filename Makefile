TARGET=fagin.loc

all:
	loc -kx locout ${TARGET}

.PHONY: run
run:
	./manifold-nexus.py run

.PHONY: clean
clean:
	rm manifold-nexus.py
	rm -rf locout archive
