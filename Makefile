all:
	cabal install \
		--bindir=./ \
		--enable-profiling

.PHONY: profile
profile:
	./fagin +RTS -p -RTS

.PHONY: test
test:
	cabal test

.PHONY: run
run:
	cabal run
