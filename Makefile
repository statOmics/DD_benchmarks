.PHONY: all benchmarks package clean clean-lupus clean-covid \
	lupus lupus-n_patients covid stagewise figures

R := Rscript --no-save --no-restore

all: benchmarks

## benchmarks: run benchmarks
benchmarks: covid lupus lupus-n_patients

lupus: package
	$(MAKE) -C benchmarks/lupus

lupus-n_patients: package
	$(MAKE) -C benchmarks/lupus-n_patients

stagewise: package
	$(MAKE) -C benchmarks/stagewise

covid: package
	$(MAKE) -C benchmarks/covid

## package: install the helper package DDCompanion
package: install.done
install.done: package/DESCRIPTION
	@echo "\nInstalling helper package"
	R CMD INSTALL package
	@touch $@

FIGURES := fig1_lupus_ncM_mock.png figS1_lupus_ncM_mock.png figS2_lupus_T4naive_mock.png \
	figS3_lupus_Bmem_mock.png figS4_lupus_mock_sce.png figS5_covid_class_switched_mock.png \
	figS6_covid_immature_mock.png figS7_covid_naive_mock.png figS9_lupus_sim_all.png \
	fig2_lupus_sim_5v5.png figS8_lupus_sim_all_alt.png figS10_covid_sim_all.png
FIGURES := $(patsubst %, figures/%, $(FIGURES))

## figures: generate figures
figures: $(FIGURES)

figures/$(FIGURES) &: figures.R
	$(R) $<

## clean: remove all generated results
clean: clean-lupus clean-covid

clean-lupus:
	$(MAKE) -C benchmarks/lupus clean

clean-covid:
	$(MAKE) -C benchmarks/covid clean
