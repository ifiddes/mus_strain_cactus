batchSystem = singleMachine
maxThreads = 10
defaultMemory = 8589934592
jobTree = .jobTree

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:${PATH}

#genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
genomes = Rattus 129S1
refGenome = C57B6J

annotationBed = /hive/groups/recon/projs/mus_strain_cactus/data/gene_check/wgEncodeGencodeBasicVM2.coding.gene-check.bed
alignmentDir = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1411/results/chained
sequenceDir = /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1411
geneCheckDir = /hive/users/markd/gencode/mus_strain_cactus/cactusMapCheck/experiments/1411/results/tracks
gencodeAttributeMap = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/data/wgEncodeGencodeAttrsVM2.attrs

all :
	cd sonLib && make
	cd jobTree && make

run : all
	if [ -d ${jobTree} ] ; then \
		rm -rf ${jobTree}; \
	fi
	python src/main.py --refGenome ${refGenome} --genomes ${genomes} --annotationBed ${annotationBed} \
	--alignmentDir ${alignmentDir} --sequenceDir ${sequenceDir}  --geneCheckDir ${geneCheckDir} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${jobTree} --overwriteDb --gencodeAttributeMap ${gencodeAttributeMap}
