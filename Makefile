all:
	cabal install

.PHONY: test
test:
	cabal test

.PHONY: run
run:
	cabal run
