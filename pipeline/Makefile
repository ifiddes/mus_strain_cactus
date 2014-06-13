.SECONDARY:
.PHONY: all clean

host=$(shell hostname)
ppid=$(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

export SHELLOPTS=pipefail
export PATH:=${PATH}:${ROOT}/extern/pycbio/bin:/hive/groups/recon/local/bin


genomes = AKRJ BALBcJ C3HHeJ C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ Rattus
refGenome = C57B6J

geneCheck = /hive/groups/recon/local/bin/gene-check
geneCheckStats = /cluster/home/markd/compbio/code/pycbio/bin/geneCheckStats
geneCheckBeds = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/tracks
sequenceDir = /hive/groups/recon/projs/mus_strain_cactus/data
alignmentDir = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-04-16.simpleChain/results/lnb_0001/chained
halFile = /hive/users/benedict/mouseHub/all.hal
filtersDir = filters
filters = ${filtersDir}/shortScaffold_noStartStop ${filtersDir}/nScan

all: $(foreach f,$(basename $(notdir ${filters})),$(foreach g,${genomes},results/$f.$g/done))

clean:
	rm -rf results/

results/%:
	mkdir -p $(dir $@)
	${filtersDir}/$(firstword $(subst ., ,$*)) \
	--refGenome ${refGenome} \
	--genome $(firstword $(subst /, , $(word 2, $(subst ., ,$*)))) \
	--geneCheckBed ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check.bed \
	--geneCheckBedDetails ${geneCheckBeds}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).coding.gene-check-details.bed \
	--alignment ${alignmentDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).chained.psl \
	--sequence ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).fa \
	--chromSizes ${sequenceDir}/$(firstword $(subst /, , $(word 2, $(subst ., ,$*)))).sizes \
	--outDir $(dir $@) > $@.${tmpExt}
	mv $@.${tmpExt} $@
