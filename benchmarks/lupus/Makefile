.PHONY: all \
	data mock_data mock_rep_data \
	results mock_results sim_results mock_results_pb sim_results_pb \
	analysis figures settings help


# Variables
R := Rscript --no-save --no-restore

RAWDATA=data-raw/GSE174188_CLUES1_adjusted.h5ad
DATA_CLEANED := data/lupus-SCE-cleaned.rds
DATA_FILTERED := data/lupus-SCE-filtered.rds
DATA := $(DATA_CLEANED) $(DATA_FILTERED)

# Subdirectories
DATADIR := data
RESDIR := results
ANALYSISDIR := analysis
FIGDIR := figures

MOCK_DATA := $(DATADIR)/lupus-SCE-mock-h5/se.rds  # 'seed' dataset for mock data

USE_CELLTYPES := B_mem ncM T4_naive # add back later
USE_ALL := $(USE_CELLTYPES) combined #TODO not sure if the combined can be removed without breaking things

MOCK_CT_SUBDIR :=$(DATADIR)/lupus-SCE-mock-ct-subsets
MOCK_CT_DATA := $(patsubst %, $(MOCK_CT_SUBDIR)/lupus-SCE-mock-%.rds, $(USE_CELLTYPES))
MOCK_REP_DATA := $(patsubst %, $(DATADIR)/lupus-SCE_list-%-mock_replicates.rds, $(USE_CELLTYPES))

PROP_DE := 0.05
_PROP_DE := $(subst .,_, $(PROP_DE))
SIM_REP_DATA := $(foreach CT,$(USE_ALL),$(foreach P,$(_PROP_DE),$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds))

PKG_STATUS := ../../install.done

# Methods used
USE_METHODS := edgeR_QP edgeR_NB \
	qbGLM_offset_squeeze bGLM \
	qbGLM qbGLM_offset

# Mock and sim results for single-cell data
#MOCK_RESULTS := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
#	$(RESDIR)/lupus-$(CT)-mock_results-sce-$(M).rds))
#SIM_RESULTS := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
#	$(foreach P,$(_PROP_DE),$(RESDIR)/lupus-$(CT)-sim_results-sce-$(M)-prop_DE-$(P).rds)))

# Mock results for pseudobulk data
MOCK_RESULTS_PB := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(RESDIR)/lupus-$(CT)-mock_results-pb-$(M).rds))
SIM_RESULTS_PB := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(foreach P,$(_PROP_DE),$(RESDIR)/lupus-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds)))

## all: generate all output files
all: analysis figures results
	@echo "\nAll jobs done."


# DATA ------------------------------------------------------------------------

## data: generate processed data
#data: $(DATA) mock_data mock_rep_data sim_rep_data

#$(DATA_CLEANED): scripts/00-anndata-to-SCE.R $(RAWDATA)
#	@echo "\nGenerating $@"
#	mkdir -p data
#	$(R) $<

#$(DATA_FILTERED): scripts/01-filter-samples.R $(DATA_CLEANED)
#	@echo "\nGenerating $@"
#	$(R) $<


## mock_data: generate mock data
#mock_data: $(MOCK_DATA) $(MOCK_CT_DATA)

#$(MOCK_DATA): $(DATA_FILTERED) scripts/02-prepare-mock.R $(PKG_STATUS)
#	@echo "\nGenerating $@"
#	mkdir -p $(DATADIR)
#	$(R) scripts/02-prepare-mock.R -v --infile $< --outfile $@

#$(MOCK_CT_DATA) &: $(MOCK_DATA) scripts/02b-prepare-mock-ct-subsets.R $(PKG_STATUS)
#	@echo "\nGenerating $(MOCK_CT_DATA)"
#	$(R) scripts/02b-prepare-mock-ct-subsets.R -v --infile $< --out_dir $(MOCK_CT_SUBDIR) \
#		--use_celltypes $(USE_CELLTYPES)

## mock_rep_data: generate mock data replicates
#mock_rep_data: $(MOCK_REP_DATA)

#$(DATADIR)/lupus-SCE_list-%-mock_replicates.rds: $(MOCK_CT_SUBDIR)/lupus-SCE-mock-%.rds \
#scripts/03-prepare-mock-replicates.R $(PKG_STATUS)
#	@echo "\nGenerating $@"
#	$(R) scripts/03-prepare-mock-replicates.R -v --infile $< --outfile $@ \
#		--n_mock_replicates 5

## sim_rep_data: generate simulated data replicates
#sim_rep_data: $(SIM_REP_DATA)

#define sim_rep_rule
#$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(subst .,_,$(P)).rds: \
#$(DATADIR)/lupus-SCE_list-$(CT)-mock_replicates.rds \
#scripts/03b-prepare-sim-replicates.R $(PKG_STATUS)
#	@echo "\nGenerating $$@"
#	$(R) scripts/03b-prepare-sim-replicates.R -v --infile $$< --outfile $$@ \
#		--prop_DE $(P)
#endef
 # --countsimQC

#$(foreach CT,$(USE_ALL),$(foreach P,$(PROP_DE),$(eval $(call sim_rep_rule))))


# RESULTS ---------------------------------------------------------------------

## results: generate results
results: mock_results_pb sim_results_pb # mock_results sim_results 

## mock_results: generate mock DE results for single-cell data
#mock_results: $(MOCK_RESULTS)

## sim_results: generate simulated DE results for single-cell data
#sim_results: $(SIM_RESULTS)

## mock_results_pb: generate mock DE results for pseudobulk data
mock_results_pb: $(MOCK_RESULTS_PB)

## sim_results_pb: generate simulated DE results for pseudobulk data
sim_results_pb: $(SIM_RESULTS_PB)

## Mock results for single-cell data
#define run_mock_rule
#$(RESDIR)/lupus-$(CT)-mock_results-sce-$(M).rds: $(DATADIR)/lupus-SCE_list-$(CT)-mock_replicates.rds \
#scripts/04-run-DE-analyses-sce.R $(PKG_STATUS)
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/04-run-DE-analyses-sce.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),$(eval $(call run_mock_rule))))

## Simulation results for single-cell data
#define run_sim_rule
#$(RESDIR)/lupus-$(CT)-sim_results-sce-$(M)-prop_DE-$(P).rds: \
#$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
#scripts/04-run-DE-analyses-sce.R $(PKG_STATUS)
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/04-run-DE-analyses-sce.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),$(foreach P,$(_PROP_DE),\
#	$(eval $(call run_sim_rule,$(CT),$(M),$(P))))))

## Mock results for pseudobulk data
define run_mock_pb_rule
$(RESDIR)/lupus-$(CT)-mock_results-pb-$(M).rds: $(DATADIR)/lupus-SCE_list-$(CT)-mock_replicates.rds \
scripts/04b-run-DE-analyses-pb.R $(PKG_STATUS)
	@echo "\nGenerating $$@ from $$<"
	mkdir -p $(RESDIR)
	$(R) scripts/04b-run-DE-analyses-pb.R -v -i $$< -o $$@ --method $(M)
endef
$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),$(eval $(call run_mock_pb_rule))))

## Simulation results for pseudobulk data
define run_sim_pb_rule
$(RESDIR)/lupus-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds: \
$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
scripts/04b-run-DE-analyses-pb.R $(PKG_STATUS)
	@echo "\nGenerating $$@ from $$<"
	mkdir -p $(RESDIR)
	$(R) scripts/04b-run-DE-analyses-pb.R -v -i $$< -o $$@ --method $(M)
endef
$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),$(foreach P,$(_PROP_DE),\
	$(eval $(call run_sim_pb_rule,$(CT),$(M),$(P))))))

# ANALYSIS --------------------------------------------------------------------

MOCK_RES_REPORT := $(ANALYSISDIR)/lupus-mock-results-sce.html
MOCK_RES_PB_REPORT := $(ANALYSISDIR)/lupus-mock-results-pb.html

## Use 1 report per celltype for the simulation results
SIM_RES_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/lupus-sim-results-$(CT)-sce.html)
SIM_RES_PB_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/lupus-sim-results-$(CT)-pb.html)

## analysis: render reports based on results
analysis: \
$(MOCK_RES_REPORT) \
$(SIM_RES_REPORTS) \
$(MOCK_RES_PB_REPORT) \
$(SIM_RES_PB_REPORTS)

## Report for single-cell mock results
$(MOCK_RES_REPORT): analysis/lupus-mock-results.Rmd $(MOCK_RESULTS) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $< \
		--out_file $@ \
		--use_celltypes $(USE_CELLTYPES) \
		--use_methods $(USE_METHODS) \
		--datatype sce

## Report for single-cell simulation results
define sim_report_rule
$(ANALYSISDIR)/lupus-sim-results-$(CT)-sce.html: analysis/lupus-sim-results.Rmd $(SIM_RESULTS) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $$< \
		--out_file $$@ \
		--use_celltypes $(CT) \
		--use_methods $(USE_METHODS) \
		--prop_DE $(PROP_DE) \
		--datatype sce
endef
$(foreach CT,$(USE_CELLTYPES),$(eval $(call sim_report_rule,$(CT))))

## Report for pseudobulk mock results
$(MOCK_RES_PB_REPORT): analysis/lupus-mock-results.Rmd $(MOCK_RESULTS_PB) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $< \
		--out_file $@ \
		--use_celltypes $(USE_CELLTYPES) \
		--use_methods $(USE_METHODS) \
		--datatype pb

## Report for pseudobulk simulation results
define sim_report_rule
$(ANALYSISDIR)/lupus-sim-results-$(CT)-pb.html: analysis/lupus-sim-results.Rmd $(SIM_RESULTS) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $$< \
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