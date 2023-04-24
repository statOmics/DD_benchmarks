.PHONY: all clean clean-data clean-mock_data clean-sim_data clean-results \
	data mock_rep_data sim_rep_data \
	results mock_results sim_results \
	analysis settings help


# Variables -------------------------------------------------------------------

R := Rscript --no-save --no-restore
PKG_STATUS := ../../install.done

# Number of patients to use for MOCK data, if missing will use all available
N_PATIENTS := 10 20 30
USE_CELLTYPES := T4_naive ncM

INPUT_DATADIR := ../lupus/data
DATADIR := data
RESDIR := results
ANALYSISDIR := analysis

PROP_DE := 0.05
_PROP_DE := $(subst .,_, $(PROP_DE))

INPUT_MOCK_DATA := \
	$(patsubst %, $(INPUT_DATADIR)/lupus-SCE_list-%-mock_replicates.rds, $(USE_CELLTYPES))
INPUT_SIM_DATA := \
	$(foreach CT,$(USE_CELLTYPES),$(foreach P,\
		$(_PROP_DE),$(INPUT_DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds))
INPUT_DATA := $(INPUT_MOCK_DATA) $(INPUT_SIM_DATA)

MOCK_REP_DATA := \
	$(foreach NP,$(N_PATIENTS),\
		$(patsubst ../lupus/data/lupus-%, data/lupus-n_patients-$(NP)-%, $(INPUT_MOCK_DATA)))
SIM_REP_DATA := \
	$(foreach NP,$(N_PATIENTS),\
		$(patsubst ../lupus/data/lupus-%, data/lupus-n_patients-$(NP)-%, $(INPUT_SIM_DATA)))
DATA := $(MOCK_REP_DATA) $(SIM_REP_DATA)


USE_METHODS := edgeR_QP edgeR_NB \
	qbGLM_offset_squeeze bGLM \
	qbGLM qbGLM_offset

# Mock and sim results for individual celltypes
MOCK_RESULTS := $(foreach NP,$(N_PATIENTS),$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(RESDIR)/lupus-n_patients-$(NP)-$(CT)-mock_results-pb-$(M).rds)))
SIM_RESULTS := $(foreach NP,$(N_PATIENTS),$(foreach M,$(USE_METHODS),$(foreach CT,$(USE_CELLTYPES),\
	$(foreach P,$(_PROP_DE),$(RESDIR)/lupus-n_patients-$(NP)-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds))))


## all: generate all output files
all: analysis results
	@echo "\nAll jobs done."


# DATA ------------------------------------------------------------------------

## data: generate processed data
#data: $(DATA)

## mock_rep_data: generate mock data replicates
#mock_rep_data: $(MOCK_REP_DATA)

#define mock_rep_rule
#$(DATADIR)/lupus-n_patients-$(NP)-SCE_list-$(CT)-mock_replicates.rds: \
#$(INPUT_DATADIR)/lupus-SCE_list-$(CT)-mock_replicates.rds \
#scripts/01-sample-patients.R
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(DATADIR)
#	$(R) scripts/01-sample-patients.R -v -i $$< -o $$@ --n_patients $(NP)
#endef
#$(foreach NP,$(N_PATIENTS),$(foreach CT,$(USE_CELLTYPES),\
#	$(eval $(call mock_rep_rule))))

## sim_rep_data: generate sim data replicates
#sim_rep_data: $(SIM_REP_DATA)

#define sim_rep_rule
#$(DATADIR)/lupus-n_patients-$(NP)-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds: \
#$(INPUT_DATADIR)/lupus-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
#scripts/01-sample-patients.R
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(DATADIR)
#	$(R) scripts/01-sample-patients.R -v -i $$< -o $$@ --n_patients $(NP)
#endef
#$(foreach NP,$(N_PATIENTS),$(foreach CT,$(USE_CELLTYPES),$(foreach P,$(_PROP_DE),\
#	$(eval $(call sim_rep_rule)))))


# RESULTS ---------------------------------------------------------------------

## results: generate results
#results: mock_results sim_results

## mock_results: generate mock DE results on individual celltypes
#mock_results: $(MOCK_RESULTS)

## sim_results: generate simulated DE results on individual celltypes
#sim_results: $(SIM_RESULTS)

## Mock results on individual celltypes
#define run_mock_rule
#$(RESDIR)/lupus-n_patients-$(NP)-$(CT)-mock_results-pb-$(M).rds: \
#$(DATADIR)/lupus-n_patients-$(NP)-SCE_list-$(CT)-mock_replicates.rds \
#scripts/02-run-DE-analyses.R
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/02-run-DE-analyses.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach NP,$(N_PATIENTS),$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),\
#	$(eval $(call run_mock_rule)))))

#define run_sim_rule
#$(RESDIR)/lupus-n_patients-$(NP)-$(CT)-sim_results-pb-$(M)-prop_DE-$(P).rds: \
#$(DATADIR)/lupus-n_patients-$(NP)-SCE_list-$(CT)-sim_replicates-prop_DE-$(P).rds \
#scripts/02-run-DE-analyses.R $(PKG_STATUS)
#	@echo "\nGenerating $$@ from $$<"
#	mkdir -p $(RESDIR)
#	$(R) scripts/02-run-DE-analyses.R -v -i $$< -o $$@ --method $(M)
#endef
#$(foreach NP,$(N_PATIENTS),$(foreach CT,$(USE_CELLTYPES),$(foreach M,$(USE_METHODS),$(foreach P,$(_PROP_DE),\
#	$(eval $(call run_sim_rule))))))


# ANALYSIS --------------------------------------------------------------------

# Use 1 report per celltype for the reports
MOCK_RES_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/lupus-n_patients-mock-results-$(CT).html)
SIM_RES_REPORTS := $(foreach CT,$(USE_CELLTYPES),$(ANALYSISDIR)/lupus-n_patients-sim-results-$(CT).html)

## analysis: render reports based on results
analysis: $(MOCK_RES_REPORTS) $(SIM_RES_REPORTS)

# Run for each celltype individually!
define mock_report_rule
$(ANALYSISDIR)/lupus-n_patients-mock-results-$(CT).html: \
analysis/lupus-n_patients-mock-results.Rmd $(MOCK_RESULTS) \
scripts/03-render-lupus-n_patients-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/03-render-lupus-n_patients-results.R -v --rmd-file $$< \
		--out_file $$@ \
		--n_patients $(N_PATIENTS) \
		--use_celltype $(CT) \
		--use_methods $(USE_METHODS)
endef
$(foreach CT,$(USE_CELLTYPES),$(eval $(call mock_report_rule,$(CT))))

define sim_report_rule
$(ANALYSISDIR)/lupus-n_patients-sim-results-$(CT).html: \
analysis/lupus-n_patients-sim-results.Rmd $(SIM_RESULTS) \
scripts/03-render-lupus-n_patients-results.R
	@echo "\nGenerating $$@"
	mkdir -p $(ANALYSISDIR)
	$(R) scripts/03-render-lupus-n_patients-results.R -v --rmd-file $$< \
		--out_file $$@ \
		--n_patients $(N_PATIENTS) \
		--use_celltype $(CT) \
		--use_methods $(USE_METHODS)
endef
$(foreach CT,$(USE_CELLTYPES),$(eval $(call sim_report_rule,$(CT))))


# Cleaning & Settings ---------------------------------------------------------

## clean-results: remove the generated results
clean-results: clean-mock_results clean-sim_results
clean-mock_results:
	rm -f $(MOCK_RESULTS)
clean-sim_results:
	rm -f $(SIM_RESULTS)

## clean: remove all generated files, except processed data
clean: clean-results
	rm -f $(MOCK_RES_REPORTS) $(SIM_RES_REPORTS)

## clean-data: remove processed data and all generated files
clean-data: clean
	rm -f $(DATA) $(MOCK_DATA)


## settings: show variables's values
settings:
	@echo R: $(R)
	@echo "\nINPUT_DATA:" $(patsubst %, "\n\t*" %, $(INPUT_DATA))
	@echo "\nDATADIR:" $(DATADIR)
	@echo "RESDIR:" $(RESDIR)
	@echo "ANALYSISDIR:" $(ANALYSISDIR)
	@echo "\nN_PATIENTS:" $(N_PATIENTS)
	@echo "USE_CELLTYPES:" $(USE_CELLTYPES)
	@echo "\nMOCK_REP_DATA:" $(patsubst %, "\n\t*" %, $(MOCK_REP_DATA))
	@echo "\nPROP_DE:" $(PROP_DE)
	@echo "SIM_REP_DATA:" $(patsubst %, "\n\t*" %, $(SIM_REP_DATA))
	@echo "\nUSE_METHODS:" $(USE_METHODS)
	@echo "\nMOCK_RESULTS:" $(patsubst %, "\n\t*" %, $(MOCK_RESULTS))
	@echo "\nSIM_RESULTS:" $(patsubst %, "\n\t*" %, $(SIM_RESULTS))
	@echo "\nMOCK_RES_REPORTS:" $(patsubst %, "\n\t*" %, $(MOCK_RES_REPORTS))
	@echo "\nSIM_RES_REPORTS:" $(patsubst %, "\n\t*" %, $(SIM_RES_REPORTS))

## help: show all commands.
help:
	@grep -h -E '^##' ${MAKEFILE_LIST} | sed -e 's/## //g' \
	| column -t -s ':'
