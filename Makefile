all:
	cabal install \
		--bindir=./ \
		--enable-profiling \
		--ghc-options="-threaded -XStrict"

.PHONY: fast
fast:
	cabal install --bindir=./ -O2

.PHONY: profile
profile:
	./fagin +RTS -p -RTS > /dev/null

.PHONY: test
test:
	cabal test

.PHONY: run
run:
	cabal run
