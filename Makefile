all:
	cabal install \
		--bindir=./ \
		--enable-profiling

.PHONY: profile
profile:
	./fagin +RTS -p -RTS > /dev/null

.PHONY: test
test:
	cabal test

.PHONY: run
run:
	cabal run
