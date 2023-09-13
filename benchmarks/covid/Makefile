.PHONY: all \
	data mock_data mock_rep_data sim_rep_data \
	results mock_results sim_results mock_results_pb sim_results_pb \
	analysis figures settings help


# Variables
R := Rscript --no-save --no-restore

DATA_CLEANED := data/covid-SCE-cleaned.rds
DATA := $(DATA_CLEANED)

#$(info $$RAWDATA is [${RAWDATA}])

# Subdirectories
DATADIR := data
RESDIR := results
ANALYSISDIR := analysis
FIGDIR := figures

MOCK_DATA := $(DATADIR)/covid-SCE-prepared.rds  # 'seed' dataset for mock data

USE_CELLTYPES := class_switched_memory_B_cell_Healthy immature_B_cell_Healthy naive_B_cell_Healthy \
	class_switched_memory_B_cell_Mild immature_B_cell_Mild naive_B_cell_Mild \
	class_switched_memory_B_cell_Moderate immature_B_cell_Moderate naive_B_cell_Moderate \
	class_switched_memory_B_cell_Critical immature_B_cell_Critical naive_B_cell_Critical

MOCK_CT_SUBDIR :=$(DATADIR)/covid-SCE-mock-ct-subsets
MOCK_CT_DATA := $(patsubst %, $(MOCK_CT_SUBDIR)/covid-SCE-mock-%.rds, $(USE_CELLTYPES))
MOCK_REP_DATA := $(patsubst %, $(DATADIR)/covid-SCE_list-%-mock_replicates.rds, $(USE_CELLTYPES))

PROP_DE := 0.05
_PROP_DE := $(subst .,_, $(PROP_DE))
SIM_REP_DATA := $(foreach CT,$(USE_CELLTYPES),$(foreach P,$(_PROP_DE),$(DATADIR)/covid-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds))

PKG_STATUS := ../../install.done

# Methods used
USE_METHODS := bGLM qbGLM \
	qbGLM_offset qbGLM_offset_squeeze \
	edgeR_NB edgeR_QP \
	edgeR_NB_optim edgeR_QP_optim

# Mock and sim  results for pseudobulk data
MOCK_RESULTS_PB := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(RESDIR)/covid-$(CT)-mock_results-pb-$(M).rds))
SIM_RESULTS_PB := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(foreach P,$(_PROP_DE),$(RESDIR)/covid-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds)))

## all: generate all output files
#all: analysis figures results
#	@echo "\nAll jobs done."


# DATA ------------------------------------------------------------------------

## data: generate processed data
#data: $(DATA) mock_data mock_rep_data sim_rep_data

#$(DATA_CLEANED): scripts/00-anndata-to-SCE.R
#	@echo "\nGenerating $@"
#	mkdir -p data
#	$(R) $<

## mock_data: generate mock data
#mock_data: $(MOCK_DATA) $(MOCK_CT_DATA)

#$(MOCK_DATA): $(DATA_CLEANED) scripts/01-prepare-mock.R $(PKG_STATUS)
#	@echo "\nGenerating $@"
#	mkdir -p $(DATADIR)
#	$(R) scripts/01-prepare-mock.R -v --infile $< --outfile $@

#$(MOCK_CT_DATA) &: $(MOCK_DATA) scripts/01b-prepare-mock-ct-subsets.R $(PKG_STATUS)
#	@echo "\nGenerating $(MOCK_CT_DATA)"
#	$(R) scripts/01b-prepare-mock-ct-subsets.R -v --infile $< --out_dir $(MOCK_CT_SUBDIR) \
#		--use_celltypes $(USE_CELLTYPES)

## mock_rep_data: generate mock data replicates
#mock_rep_data: $(MOCK_REP_DATA)

#$(DATADIR)/covid-SCE_list-%-mock_replicates.rds: $(MOCK_CT_SUBDIR)/covid-SCE-mock-%.rds \
#scripts/02-prepare-mock-replicates.R $(PKG_STATUS)
#	@echo "\nGenerating $@"
#	$(R) scripts/02-prepare-mock-replicates.R -v --infile $< --outfile $@ \
#		--n_mock_replicates 5

## sim_rep_data: generate simulated data replicates
#sim_rep_data: $(SIM_REP_DATA)

#define sim_rep_rule
#$(DATADIR)/covid-SCE_list-$(CT)-sim_replicates-prop_DE-$(subst .,_,$(P)).rds: \
#$(DATADIR)/covid-SCE_list-$(CT)-mock_replicates.rds \
#scripts/02b-prepare-sim-replicates.R $(PKG_STATUS)
#	@echo "\nGenerating $$@"
#	$(R) scripts/02b-prepare-sim-replicates.R -v --infile $$< --outfile $$@ \
#		--prop_DE $(P)
#endef

#$(foreach CT,$(USE_CELLTYPES),$(foreach P,$(PROP_DE),$(eval $(call sim_rep_rule))))


# RESULTS ---------------------------------------------------------------------

## results: generate results
#results: mock_results_pb sim_results_pb

## mock_results_pb: generate mock DE results for pseudobulk data
#mock_results_pb: $(MOCK_RESULTS_PB)

## sim_results_pb: generate simulated DE results for pseudobulk data
#sim_results_pb: $(SIM_RESULTS_PB)

## Mock results for pseudobulk data
#define run_mock_pb_rule
#$(RESDIR)/covid-$(CT)-mock_results-pb-$(M).rds: $(DATADIR)/covid-SCE_list-$(CT)-mock_replicates.rds \
#scripts/03-run-DE-analyses-pb.R $(PKG_STATUS)
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/03-run-DE-analyses-pb.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),$(eval $(call run_mock_pb_rule))))

## Simulation results for pseudobulk data
#define run_sim_pb_rule
#$(RESDIR)/covid-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds: $(DATADIR)/covid-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
#scripts/03-run-DE-analyses-pb.R $(PKG_STATUS)
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/03-run-DE-analyses-pb.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),$(foreach P,$(_PROP_DE),\
#	$(eval $(call run_sim_pb_rule,$(CT),$(M),$(P))))))

# ANALYSIS --------------------------------------------------------------------

MOCK_RES_PB_REPORT := $(ANALYSISDIR)/covid-mock-results-pb.html

## Use 1 report per celltype for the simulation results
SIM_RES_PB_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/covid-sim-results-$(CT)-pb.html)

## analysis: render reports based on results
analysis: \
$(MOCK_RES_PB_REPORT) \
$(SIM_RES_PB_REPORTS)

## Report for pseudobulk mock results
$(MOCK_RES_PB_REPORT): analysis/covid-mock-results.Rmd $(MOCK_RESULTS_PB) \
scripts/04-render-covid-results.R
	@echo "\nGenerating $@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/04-render-covid-results.R -v --rmd-file $< \
		--out_file $@ \
		--use_celltypes $(USE_CELLTYPES) \
		--use_methods $(USE_METHODS) \
		--datatype pb

## Report for pseudobulk simulation results
define sim_report_rule
$(ANALYSISDIR)/covid-sim-results-$(CT)-pb.html: analysis/covid-sim-results.Rmd $(SIM_RESULTS) \
scripts/04-render-covid-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/04-render-covid-results.R -v --rmd-file $$< \
		--out_file $$@ \
		--use_celltypes $(CT) \
		--use_methods $(USE_METHODS) \
		--prop_DE $(PROP_DE) \
		--datatype pb
endef
$(foreach CT,$(USE_CELLTYPES),$(eval $(call sim_report_rule,$(CT))))

## settings: show variables's values
settings:
	@echo R: $(R)
	@echo "\nRAWDATA (input data):" $(RAWDATA)
	@echo "DATA (processed data):" $(DATA)
	@echo "\nDATADIR:" $(DATADIR)
	@echo "RESDIR:" $(RESDIR)
	@echo "ANALYSISDIR:" $(ANALYSISDIR)
	@echo "FIGDIR:" $(FIGDIR)
	@echo "\nUSE_CELLTYPES:" $(USE_CELLTYPES)
	@echo "MOCK_DATA:" $(MOCK_DATA)
	@echo "MOCK_CT_DATA:" $(patsubst %, "\n\t*" %, $(MOCK_CT_DATA))
	@echo "MOCK_REP_DATA:" $(patsubst %, "\n\t*" %, $(MOCK_REP_DATA))
	@echo "\nPROP_DE:" $(PROP_DE)
	@echo "SIM_REP_DATA:" $(patsubst %, "\n\t*" %, $(SIM_REP_DATA))
	@echo "\nUSE_METHODS:" $(USE_METHODS)
	@echo "\nMOCK_RESULTS:" $(patsubst %, "\n\t*" %, $(MOCK_RESULTS))
	@echo "\nSIM_RESULTS:" $(patsubst %, "\n\t*" %, $(SIM_RESULTS))
	@echo "\nMOCK_RES_REPORT:" $(MOCK_RES_REPORT)
	@echo "\nMOCK_RES_PB_REPORT:" $(MOCK_RES_PB_REPORT)
	@echo "SIM_RES_REPORTS:" $(patsubst %, "\n\t*" %, $(SIM_RES_REPORTS))
	@echo "\nFIGURES:" $(patsubst %, "\n\t*" %, $(FIGURES))

## help: show all commands.
help:
	@grep -h -E '^##' ${MAKEFILE_LIST} | sed -e 's/## //g' \
	| column -t -s ':'
