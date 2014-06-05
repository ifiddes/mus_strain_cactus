.SECONDARY:
.PHONY: all clean

host=$(shell hostname)
ppid=$(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

export SHELLOPTS=pipefail
export PATH:=${PATH}:${ROOT}/extern/pycbio/bin:/hive/groups/recon/local/bin


genomes = AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ Rattus

geneCheck = /hive/groups/recon/local/bin/gene-check
geneCheckStats = /cluster/home/markd/compbio/code/pycbio/bin/geneCheckStats
geneCheckResults = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/geneCheck/
genomeDir = /hive/groups/recon/projs/mus_strain_cactus/data
filtersDir = filters
filters = $(shell find ${filtersDir}/ -perm 755 -type f)


all: $(foreach f,$(basename $(notdir ${filters})),$(foreach g,${genomes},results/$f.$g/stdout))

clean:
	rm -rf results/

results/%:
	mkdir -p $$(dirname $@)
	echo $@
	${filtersDir}/$(firstword $(subst ., ,$*)) gene_check_bed gene_check_bed_details genome ${genomeDir} > $@.${tmpExt}
	mv $@.${tmpExt} $@
