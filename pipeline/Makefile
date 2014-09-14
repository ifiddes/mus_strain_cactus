.SECONDARY:
.PHONY: all clean check-release

host=$(shell hostname)
ppid=$(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

export SHELLOPTS=pipefail
export PATH:=${PATH}:${ROOT}/extern/pycbio/bin:/hive/groups/recon/local/bin

genomes_1302 = Rattus AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ
genomes_1405 = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

dataDir = /hive/groups/recon/projs/mus_strain_cactus/data
geneCheck = /hive/groups/recon/local/bin/gene-check
geneCheckStats = /cluster/home/markd/compbio/code/pycbio/bin/geneCheckStats
originalGeneCheckBeds = ${dataDir}/gene_check
geneCheckBeds_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/tracks
geneCheckBeds_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/tracks
augustus_beds= /hive/groups/recon/projs/mus_strain_cactus/experiments/016.augustus_to_genepred
geneCheckBeds_c6tr_cgp = ${augustus_beds}/c6tr_cgp
geneCheckBeds_c6tr_mea = ${augustus_beds}/c6tr_mea
geneCheckBeds_chr6_cgp = ${augustus_beds}/chr6_cgp
geneCheckBeds_chr6_mea = ${augustus_beds}/chr6_mea
alignmentDir_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/chained
alignmentDir_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/chained
alignmentDir_augustus = ${augustus_beds}/fake_alignments
halFile = /hive/users/benedict/mouseHub/all.hal
filtersDir = filters
stats = src/geneCheckStatCreator.py
filters = ${filtersDir}/metaFilter

ifeq ($(release),1302)
	genomes = $(genomes_1302)
	sequenceDir = ${dataDir}/assembly_rel_${release}
	geneCheckBeds = $(geneCheckBeds_1302)
	alignmentDir = $(alignmentDir_1302)
	originalGeneCheckBed = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check.bed
	originalGeneCheckBedDetails = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed
	metaFilterMode = transmap
else ifeq ($(release),1405)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_${release}
	geneCheckBeds = $(geneCheckBeds_1405)
	alignmentDir = $(alignmentDir_1405)
	originalGeneCheckBed = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check.bed
	originalGeneCheckBedDetails = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed
	metaFilterMode = transmap
else ifeq ($(release),c6tr_cgp)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_c6tr_cgp)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
	metaFilterMode = augustus
else ifeq ($(release),c6tr_mea)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_c6tr_mea)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
	metaFilterMode = augustus
else ifeq ($(release),chr6_cgp)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_chr6_cgp)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
	metaFilterMode = augustus
else ifeq ($(release),chr6_mea)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_chr6_mea)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
	metaFilterMode = augustus
else
	$(error release is unrecognized. specify 1302, 1405 or some augustus release)
endif

all: check-release $(foreach f,$(basename $(notdir ${filters})),$(foreach g,${genomes},results_${release}/$f.$g/done))

check-release:
ifndef release
	$(error release is undefined. specify 1302, 1405 or some augustus release)
endif

clean:
	echo 'youll need to run your own clean.'

results_$(release)/%:
	mkdir -p $(dir $@)
	${filtersDir}/$(firstword $(subst ., ,$*)) \
	--refGenome ${refGenome} \
	--genome $(firstword $(subst /, , $(word 2, $(subst ., ,$*)))) \
	--geneCheckBed ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check.bed \
	--geneCheckBedDetails ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check-details.bed \
	--originalGeneCheckBed ${originalGeneCheckBed} \
	--originalGeneCheckBedDetails ${originalGeneCheckBedDetails} \
	--alignment ${alignmentDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).chained.psl \
	--sequence ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).fa \
	--refSequence ${sequenceDir}/C57B6J.fa \
	--chromSizes ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).sizes \
  --mode ${metaFilterMode} \
	--outDir $(dir $@) &> $@.${tmpExt}
	${stats} --geneCheck $(dir $@)out.bed \
	--geneCheckDetails $(dir $@)out_details.bed \
	--out $(dir $@)stats.xml &> $(dir $@)stats.out
	mv $@.${tmpExt} $@
