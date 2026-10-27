# Convenience targets. Every target is a thin wrapper over scripts/.
.PHONY: help setup smoke figures tables check-release track-a track-a2 track-b track-c

help:
	@echo "make setup          bootstrap environments and fetch public inputs"
	@echo "make smoke          tests that need no data and no GPU"
	@echo "make figures        regenerate docs/figures (method schematics)"
	@echo "make tables         regenerate docs/tables"
	@echo "make check-release  audit the repository before a push"
	@echo "make track-b        ADAR guide design (plan only)"
	@echo "make track-a LIB=/path/to/enamine_real"

setup:
	scripts/bootstrap.sh all
	scripts/fetch_public_inputs.sh all

smoke:
	scripts/run_all.sh smoke

figures:
	python docs/figures/make_figures.py

tables:
	python docs/make_tables.py

check-release:
	scripts/check_release.sh

track-a:
	scripts/run_track_a.sh screen --library "$(LIB)"

track-a2:
	scripts/run_track_a2.sh

track-b:
	scripts/run_track_b.sh preflight
	scripts/run_track_b.sh plan

track-c:
	scripts/run_track_c.sh prepare
