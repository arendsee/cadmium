all:
	cabal install \
		--bindir=./ \
		--enable-profiling \
		--enable-benchmarks

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

.PHONY: bench
bench:
	cabal bench --benchmark-options="-o report.html" +RTS -T

.PHONY: clean
clean:
	rm -f *.aux *.ps *.hp *.prof
