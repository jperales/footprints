
import os
import os.path
import re
import glob
import yaml

YAML = glob.glob("index/*/*.yaml")
ACCESSION = [os.path.splitext(os.path.basename(k))[0] for k in YAML]

rule all:
	input:
		#expand("data/normalized/{acc}.RData", acc=ACCESSION)
		#rdat = "data/expr.RData"
		"data/zscores.RData",
		"data/zscores_commonNULL.RData"


rule normalize:
	input:
		script = "workflow/scripts/data_normalization.R",
		#NOTE: Although not used, we state the input file for file stamp
		#TIP: get by substring,  yaml = lambda wildcards: [k for k in YAML if wildcards.acc in k][0]
		yaml = lambda wildcards: YAML[ACCESSION.index(wildcards.acc)]
	output:
		rdat = "data/normalized/{acc}.RData"
	params:
		acc = "{acc}"
	message:
		"Normalizing Accessions"
	shell:
		'set +eu '
		' && (test -d data/normalized/ || mkdir data/normalized) '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/data '
		" && $CONDA_PREFIX/bin/Rscript {input.script}"
		# Inputs
		" {params.acc} {output.rdat}"

rule merge:
	input: 
		script = "workflow/scripts/data_expr.R",
		rdats = expand("data/normalized/{acc}.RData", acc=ACCESSION)
	output:
		rdat = "data/expr.RData"
	message:
		"Merging normalized datasets"
	shell:
		'set +eu '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/data '
		" && $CONDA_PREFIX/bin/Rscript {input.script}"
		# Inputs
		" {input.rdats}"

rule zscore:
	input:
		script = "workflow/scripts/data_zscores.R",
		rdat = "data/expr.RData"
	output:
		rdat = "data/zscores.RData"
	message:
		"Z-scoring expression"
	shell:
		'set +eu '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/data '
		" && $CONDA_PREFIX/bin/Rscript {input.script} {input.rdat} {output.rdat}"

rule zscoreNULL:
	input:
		script = "workflow/scripts/data_zscores_commonNULL.R",
		rdat = "data/expr.RData"
	output:
		rdat = "data/zscores_commonNULL.RData"
	message:
		"Z-scoring expression"
	shell:
		'set +eu '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/data '
		" && $CONDA_PREFIX/bin/Rscript {input.script} {input.rdat} {output.rdat}"

