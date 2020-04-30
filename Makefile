.PHONY: simulation matching skymap

all: simulation matching skymap

simulation:
	$(MAKE) -C $@

matching:
	$(MAKE) -C $@

skymap:
	cd $@ && ./makeskymap.py
