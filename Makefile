.PHONY: data empirics simulations all

data:
	Rscript scripts/run_data.R

empirics:
	Rscript scripts/run_empirics.R

simulations:
	Rscript scripts/run_simulations.R

all:
	Rscript scripts/run_all.R
