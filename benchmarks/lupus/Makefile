.PHONY: all \
	data mock_data mock_rep_data \
	results mock_results sim_results \
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

USE_CELLTYPES := B_mem ncM
USE_ALL := $(USE_CELLTYPES) combined #TODO not sure if the combined can be removed without breaking things

MOCK_CT_SUBDIR :=$(DATADIR)/lupus-SCE-mock-ct-subsets
MOCK_CT_DATA := $(patsubst %, $(MOCK_CT_SUBDIR)/lupus-SCE-mock-%.rds, $(USE_CELLTYPES))
MOCK_REP_DATA := $(patsubst %, $(DATADIR)/lupus-SCE_list-%-mock_replicates.rds, $(USE_CELLTYPES))

PROP_DE := 0.05
_PROP_DE := $(subst .,_, $(PROP_DE))
SIM_REP_DATA := $(foreach CT,$(USE_ALL),$(foreach P,$(_PROP_DE),$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds))

PKG_STATUS := ../../install.done

# Methods used for both type of analyses # this is an artifact from Milan's code, but I keep it because we can distinguish between sc and pb later
USE_METHODS_BOTH := muscat

# Methods used for the "single celltype" analyses only
USE_METHODS := $(USE_METHODS_BOTH) \
	qbGLM_offset_squeeze bGLM \
	qbGLM qbGLM_offset

# Mock and sim results for individual celltypes
MOCK_RESULTS := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(RESDIR)/lupus-$(CT)-mock_results-$(M).rds))
SIM_RESULTS := $(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(foreach P,$(_PROP_DE),$(RESDIR)/lupus-$(CT)-sim_results-$(M)-prop_DE-$(P).rds)))

## all: generate all output files
all: analysis figures results
	@echo "\nAll jobs done."


# DATA ------------------------------------------------------------------------

## data: generate processed data
data: $(DATA) mock_data mock_rep_data sim_rep_data

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
results: mock_results sim_results

## mock_results: generate mock DE results on individual celltypes
mock_results: $(MOCK_RESULTS)

## sim_results: generate simulated DE results on individual celltypes
sim_results: $(SIM_RESULTS)

# Mock results on individual celltypes
define run_mock_rule
$(RESDIR)/lupus-$(CT)-mock_results-$(M).rds: $(DATADIR)/lupus-SCE_list-$(CT)-mock_replicates.rds \
scripts/04-run-DE-analyses.R $(PKG_STATUS)
	@echo "\nGenerating $$@ from $$<"
	mkdir -p $(RESDIR)
	$(R) scripts/04-run-DE-analyses.R -v -i $$< -o $$@ --method $(M)
endef
$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),$(eval $(call run_mock_rule))))

define run_sim_rule
$(RESDIR)/lupus-$(CT)-sim_results-$(M)-prop_DE-$(P).rds: \
$(DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
scripts/04-run-DE-analyses.R $(PKG_STATUS)
	@echo "\nGenerating $$@ from $$<"
	mkdir -p $(RESDIR)
	$(R) scripts/04-run-DE-analyses.R -v -i $$< -o $$@ --method $(M)
endef
$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),$(foreach P,$(_PROP_DE),\
	$(eval $(call run_sim_rule,$(CT),$(M),$(P))))))

# ANALYSIS --------------------------------------------------------------------

MOCK_RES_REPORT := $(ANALYSISDIR)/lupus-mock-results.html

# Use 1 report per celltype for the simulation results
SIM_RES_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/lupus-sim-results-$(CT).html)

## analysis: render reports based on results
analysis: \
$(MOCK_RES_REPORT) \
$(SIM_RES_REPORTS)

$(MOCK_RES_REPORT): analysis/lupus-mock-results.Rmd $(MOCK_RESULTS) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $< \
		--out_file $@ \
		--use_celltypes $(USE_CELLTYPES) \
		--use_methods $(USE_METHODS)

# Run for each celltype individually!
define sim_report_rule
$(ANALYSISDIR)/lupus-sim-results-$(CT).html: analysis/lupus-sim-results.Rmd $(SIM_RESULTS) \
scripts/05-render-lupus-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/05-render-lupus-results.R -v --rmd-file $$< \
		--out_file $$@ \
		--use_celltypes $(CT) \
		--use_methods $(USE_METHODS) \
		--prop_DE $(PROP_DE)
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
	@echo "SIM_RES_REPORTS:" $(patsubst %, "\n\t*" %, $(SIM_RES_REPORTS))
	@echo "\nFIGURES:" $(patsubst %, "\n\t*" %, $(FIGURES))

## help: show all commands.
help:
	@grep -h -E '^##' ${MAKEFILE_LIST} | sed -e 's/## //g' \
	| column -t -s ':'
