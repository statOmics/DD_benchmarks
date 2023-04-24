.PHONY: all lupus ledergor benchmarks package clean clean-lupus clean-ledergor \
	lupus lupus-n_patients ledergor figures

R := Rscript --no-save --no-restore

all: benchmarks

## benchmarks: run benchmarks
benchmarks: ledergor lupus lupus-n_patients

lupus: package
	$(MAKE) -C benchmarks/lupus

lupus-n_patients: package
	$(MAKE) -C benchmarks/lupus-n_patients

ledergor: package
	$(MAKE) -C benchmarks/ledergor

## package: install the helper package DDCompanion
package: install.done
install.done: package/DESCRIPTION
	@echo "\nInstalling helper package"
	R CMD INSTALL package
	@touch $@

FIGURES := figure1.pdf figure2.pdf figure3.pdf figure4.pdf figure5.pdf
FIGURES := $(patsubst %, figures/%, $(FIGURES))

## figures: generate figures
figures: $(FIGURES)

figures/$(FIGURES) &: figures.R
	$(R) $<

## clean: remove all generated results
clean: clean-lupus clean-ledergor

clean-lupus:
	$(MAKE) -C benchmarks/lupus clean

clean-ledergor:
	$(MAKE) -C benchmarks/ledergor clean
