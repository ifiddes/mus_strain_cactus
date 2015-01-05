.SECONDARY:
.PHONY: all clean check-release

host=$(shell hostname)
ppid=$(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

export SHELLOPTS=pipefail
export PATH=:./sonLib/bin:./submodules/jobTree/bin:${PATH}
export PYTHONPATH=:./:./submodules:${PYTHONPATH}

genomes_1302 = Rattus AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ
genomes_1405 = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
genomes_1409 = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
genomes_1411 = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

originalGeneCheckBeds = ${dataDir}/gene_check
geneCheckBeds_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/tracks
geneCheckBeds_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/tracks
geneCheckBeds_1409 = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1409/results/tracks
geneCheckBeds_1411 = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1411/results/tracks
augustus_beds= /hive/groups/recon/projs/mus_strain_cactus/experiments/016.augustus_to_genepred
geneCheckBeds_c6tr_cgp = ${augustus_beds}/c6tr_cgp
geneCheckBeds_c6tr_mea = ${augustus_beds}/c6tr_mea
geneCheckBeds_chr6_cgp = ${augustus_beds}/chr6_cgp
geneCheckBeds_chr6_mea = ${augustus_beds}/chr6_mea
alignmentDir_1302 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/chained
alignmentDir_1405 = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/results/lnb_0001/chained
alignmentDir_1409 = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1409/results/chained
alignmentDir_1411 = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1411/results/chained

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
else ifeq ($(release),1409)
	genomes = $(genomes_1409)
	sequenceDir = ${dataDir}/assembly_rel_${release}
	geneCheckBeds = $(geneCheckBeds_1409)
	alignmentDir = $(alignmentDir_1409)
	originalGeneCheckBed = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check.bed
	originalGeneCheckBedDetails = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed
else ifeq ($(release),1411)
	genomes = $(genomes_1411)
	sequenceDir = ${dataDir}/assembly_rel_${release}
	geneCheckBeds = $(geneCheckBeds_1411)
	alignmentDir = $(alignmentDir_1411)
	originalGeneCheckBed = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check.bed
	originalGeneCheckBedDetails = ${originalGeneCheckBeds}/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed
else ifeq ($(release),c6tr_cgp)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_c6tr_cgp)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
else ifeq ($(release),c6tr_mea)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_c6tr_mea)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
else ifeq ($(release),chr6_cgp)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_chr6_cgp)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
else ifeq ($(release),chr6_mea)
	genomes = $(genomes_1405)
	sequenceDir = ${dataDir}/assembly_rel_1405
	geneCheckBeds = $(geneCheckBeds_chr6_mea)
	alignmentDir = $(alignmentDir_augustus)
	originalGeneCheckBed = ${augustus_beds}/originalGeneCheckBed.bed
	originalGeneCheckBedDetails = ${augustus_beds}/originalGeneCheckBedDetails.bed
else
	$(error release is unrecognized. specify 1302, 1405, 1409, 1411, or some augustus release)
endif

all: check-release $(foreach f,$(basename $(notdir ${filters})),$(foreach g,${genomes},results_${release}/$f.$g/done))

check-release:
ifndef release
	$(error release is undefined. specify 1302, 1405 or some augustus release)
endif

clean:
	echo 'youll need to run your own clean.'

results_$(release)/%:
	python src/main.py --refGenome ${refGenome} --genomes ${genomes} --geneCheckDir ${geneCheckBeds} \
	--alignmentDir ${alignmentDir} --sequenceDir ${sequenceDir} --originalGeneCheckBed ${originalGeneCheckBed} \
	--originalGeneCheckBedDetails ${originalGeneCheckBedDetails} --outDb results_${release}.db
