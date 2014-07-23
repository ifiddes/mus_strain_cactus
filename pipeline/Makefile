.SECONDARY:
.PHONY: all clean check-release

host=$(shell hostname)
ppid=$(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

export SHELLOPTS=pipefail
export PATH:=${PATH}:${ROOT}/extern/pycbio/bin:/hive/groups/recon/local/bin

genomes_1302 = AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ Rattus
genomes_1405 = AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ Rattus 129S1 AJ JASTEiJ LPJ NODSHiLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

dataDir = /hive/groups/recon/projs/mus_strain_cactus/data
geneCheck = /hive/groups/recon/local/bin/gene-check
geneCheckStats = /cluster/home/markd/compbio/code/pycbio/bin/geneCheckStats
originalGeneCheckBeds = ${dataDir}/gene_check
geneCheckBeds_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/tracks
geneCheckBeds_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/tracks
sequenceDir = ${dataDir}/assembly_rel_${release}
alignmentDir_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/chained
alignmentDir_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/chained
halFile = /hive/users/benedict/mouseHub/all.hal
filtersDir = filters
stats = src/geneCheckStatCreator.py
filters = ${filtersDir}/metaFilter

ifeq ($(release),1302)
	genomes = $(genomes_1302)
	geneCheckBeds = $(geneCheckBeds_1302)
	alignmentDir = $(alignmentDir_1302)
else
	genomes = $(genomes_1405)
	geneCheckBeds = $(geneCheckBeds_1405)
	alignmentDir = $(alignmentDir_1405)
endif

all: check-release $(foreach f,$(basename $(notdir ${filters})),$(foreach g,${genomes},results_${release}/$f.$g/done))

check-release:
ifndef release
	$(error release is undefined. specify 1302 or 1405)
endif

clean:
	rm -rf results*/

results_$(release)/%: check-release
	mkdir -p $(dir $@)
	${filtersDir}/$(firstword $(subst ., ,$*)) \
	--refGenome ${refGenome} \
	--genome $(firstword $(subst /, , $(word 2, $(subst ., ,$*)))) \
	--geneCheckBed ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check.bed \
	--geneCheckBedDetails ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check-details.bed \
	--originalGeneCheckBed ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check.bed \
	--originalGeneCheckBedDetails ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed \
	--alignment ${alignmentDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).chained.psl \
	--sequence ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).fa \
	--chromSizes ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).sizes \
	--outDir $(dir $@) &> $@.${tmpExt}
	${stats} --geneCheck $(dir $@)out.bed \
	--geneCheckDetails $(dir $@)out_details.bed \
	--out $(dir $@)stats.xml &> $(dir $@)stats.out
	mv $@.${tmpExt} $@
