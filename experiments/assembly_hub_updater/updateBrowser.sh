#!/bin/bash -e
# update the mouse browser.
# not robust in the slightest.

if [ -z $1 ]; then
  echo "Error, you must specify the release number as the first argument, either 1302 or 1405"
  exit 1
fi
release=$1

if [ "$release" != "1302" ] && [ "$release" != "1405" ]; then
  echo "Release must be either 1302 or 1405"
fi

# Rattus is put at the front of the list (out of alphabetical order) because
# it takes an extremely long time to run. I'd rather get it done first than
# be sitting and waiting for it to finish at the very end.
if [ "$release" == "1302" ]; then
  genomes="Rattus Anc0 Anc1 Anc2 Anc3 Anc4 Anc5 Anc6 Anc7 Anc8 AKRJ BALBcJ C3HHeJ C57B6J C57B6NJ CBAJ DBA2J FVBNJ NZOHlLtJ"
else
  genomes="Rattus Anc00 Anc01 Anc02 Anc03 Anc04 Anc05 Anc06 Anc07 Anc08 Anc09 Anc10 Anc11 Anc12 Anc13 Anc14 Anc15 Anc16 129S1 AJ AKRJ BALBcJ CASTEiJ C3HHeJ C57B6J C57B6NJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ"
fi

# python
PYTHONPATH=$PYTHONPATH:/hive/users/dearl/msca/proj/src/progressiveCactus/submodules
PATH=$PATH:/hive/users/dearl/msca/proj/src/progressiveCactus/submodules/hal/bin

echo "# Creating a hub for release $release"
mkdir -p ../mouseBrowser_$release

echo '# Updating bed files'
mkdir -p ./bedDirs_$release
/hive/users/dearl/msca/mus_strain_cactus/pipeline/src/result2bedDirs.py --outDir ./bedDirs_$release /hive/users/dearl/msca/mus_strain_cactus/pipeline/results_$release/metaFilter.*

echo '# Bringing over original input gene-check beds'
mkdir -p bedDirs_$release/input_geneCheck/C57B6J
mkdir -p bedDirs_$release/input_geneCheck_details/C57B6J
cp /hive/groups/recon/projs/mus_strain_cactus/data/gene_check/wgEncodeGencodeBasicVM2.coding.gene-check.bed /hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/input_geneCheck/C57B6J/C57B6J.wgEncodeGencodeBasicVM2.coding.gene-check.bed
cp /hive/groups/recon/projs/mus_strain_cactus/data/gene_check/wgEncodeGencodeBasicVM2.coding.gene-check-details.bed /hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/input_geneCheck_details/C57B6J/C57B6J.wgEncodeGencodeBasicVM2.coding.gene-check-details.bed

echo '# Removing old symlinks'
for d in refGene knownGene wgEncodeGencodeCompVM2 wgEncodeGencodeCompVM2_CDS; do
  rm -f /cluster/home/dearl/msca/myMouseBrowser/browser_$release/liftoverbed/$d/$d
done

if [ ! -d "/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/refGene" ]; then
  echo "# Extracting refGene, knownGene, wgEncodeGencodeComp data."
  mkdir -p liftover_$release
  hgsql -e "select * from refGene" mm10 | cut -f 2- > ./liftover_$release/refGene.gp &
  hgsql -e "select * from knownGene" mm10  | cut -f 2- > ./liftover_$release/knownGene.gp &
  hgsql -e "select * from wgEncodeGencodeCompVM2" mm10  | cut -f 2- > ./liftover_$release/wgEncodeGencodeCompVM2.gp &
  wait
  python ./ngan_scripts/gp2bed.py liftover_$release/refGene.gp liftover_$release/refGene.bed &
  python ./ngan_scripts/knowngene2bed.py liftover_$release/knownGene.gp liftover_$release/knownGene.bed &
  python ./ngan_scripts/gp2bed.py liftover_$release/wgEncodeGencodeCompVM2.gp liftover_$release/wgEncodeGencodeCompVM2.bed &
  wait
  echo "# Performing liftover (this takes several hours)"
  for g in $genomes; do
    for b in refGene knownGene wgEncodeGencodeCompVM2; do
      ##########
      # This block holds the background process count in check
      while [ $(jobs | wc -l) -gt 8 ]; do
        sleep 5
      done
      /hive/groups/recon/projs/mus_strain_cactus/src/progressiveCactus/submodules/hal/bin/halLiftover /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/msca.hal C57B6J liftover_$release/$b.bed $g liftover_$release/$b.$g.bed_tmp --keepExtra --outBedVersion 13 --tab &> liftover_$release/_log.liftover.$b.$g &
    done
  done
  wait
  echo "# Sorting lift over beds"
  for g in $genomes; do
    for b in refGene knownGene wgEncodeGencodeCompVM2; do
      bedSort liftover_$release/$b.$g.bed_tmp liftover_$release/$b.$g.bed.srt &
    done
    wait
    for b in refGene knownGene wgEncodeGencodeCompVM2; do
      rm liftover_$release/$b.$g.bed_tmp
      mv liftover_$release/$b.$g.bed.srt liftover_$release/$b.$g.bed
    done
  done
  echo "# Creating big beds from liftover"
  for b in refGene knownGene wgEncodeGencodeCompVM2; do
    mkdir -p /hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J
    cp asfiles/$b.C57B6J.as /hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J/C57B6J.as
  done
  for g in $genomes; do
    for b in refGene knownGene wgEncodeGencodeCompVM2; do
      # these will lack thick-thin distinction
      bedToBigBed -tab -type=bed12+1 -extraIndex=name,commonName -as=/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J/C57B6J.as liftover_$release/$b.$g.bed /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/$g.sizes /hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J/$g.bb &> _log.bedToBigBed.$b.$g  &
    done
    wait
  done
  echo "# Placing the original C57B6J big beds into the browser"
  for b in refGene knownGene wgEncodeGencodeCompVM2; do
    # we do this to get the thick-thin distinction
    bedSort liftover_$release/$b.bed liftover_$release/$b.bed.srt
    bedToBigBed -tab -type=bed12+1 -extraIndex=name,commonName -as=/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J/C57B6J.as liftover_$release/$b.bed.srt /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/C57B6J.sizes /hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/$b/C57B6J/C57B6J.bb &> _log.bedToBigBed.$b.C57B6J  &
  done
  wait
fi

if [ ! -f /hive/users/dearl/msca/mouseBrowser_$release/lod.txt ]; then
  echo '# Generating hal LOD files (this can take a day or so)'
  pushd /hive/users/dearl/msca/proj/src/progressiveCactus/submodules/ > /dev/null
  OLD_PATH=$PATH
  PATH=$PATH:/hive/users/dearl/msca/proj/src/progressiveCactus/submodules/hal/bin
  /hive/users/dearl/msca/proj/src/progressiveCactus/submodules/hal/bin/halLodInterpolate.py /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/msca.hal /hive/users/dearl/msca/mouseBrowser_$release/lod.txt.tmp --outHalDir /hive/users/dearl/msca/mouseBrowser_$release/lod --numProc 12
  mv /hive/users/dearl/msca/mouseBrowser_$release/lod.txt.tmp /hive/users/dearl/msca/mouseBrowser_$release/lod.txt
  PATH=$OLD_PATH
  popd > /dev/null
fi

echo '# Removing old job tree'
rm -rf jt_assembly_hub_$release

echo '# Launching hal2assemblyHub.py'
# create an assembly hub
/hive/users/dearl/msca/proj/src/progressiveCactus/submodules/hal/bin/hal2assemblyHub.py /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/msca.hal \
/hive/users/dearl/msca/myMouseBrowser/browser_$release \
--hub=msca_$release --shortLabel=MouseStrain_$release \
--longLabel="Mouse Strain Comparative Assembly Hub release $release" --email='benedict@soe.ucsc.edu' \
--url="http://hgwdev.sdsc.edu/~benedict/mouseBrowser_$release" \
--twobitdir=/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release \
--lod --lodTxtFile=/hive/users/dearl/msca/mouseBrowser_$release/lod.txt \
--lodDir=/hive/users/dearl/msca/mouseBrowser_$release/lod \
--finalBigBedDirs=/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/refGene,/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/knownGene,/hive/users/dearl/msca/myMouseBrowser/bigBedDirs_$release/wgEncodeGencodeCompVM2 \
--bedDirs=/hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/metaFilter,/hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/metaFilter_details,/hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/input_geneCheck,/hive/users/dearl/msca/myMouseBrowser/bedDirs_$release/input_geneCheck_details \
--tabBed \
--jobTree=./jt_assembly_hub_$release \
--batchSystem=singleMachine \
--stats \
--maxThreads=16 \
--logInfo \
--noBedLiftover

echo '# Correcting genomes.txt'
perl -ple "s/orderKey 4800/orderKey 4800\ndescription rel_$release/" -i browser_$release/genomes.txt

echo '# Copying genome files.'
# copy over genome files
for g in $genomes documentation; do
  mkdir -p ../mouseBrowser_$release/$g/
  cp -r browser_$release/$g/* ../mouseBrowser_$release/$g/ &
done
wait

if [ ! -f ../mouseBrowser_$release/msca.hal ]; then
  echo '# Copying main hal file.'
  cp /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/msca.hal ../mouseBrowser_$release/msca.hal
fi

echo '# Copying bed files.'
# copy over bed files
for d in metaFilter metaFilter_details input_geneCheck input_geneCheck_details; do
  mkdir -p ../mouseBrowser_$release/liftoverbed/$d;
  cp -r browser_$release/liftoverbed/$d/* ../mouseBrowser_$release/liftoverbed/$d/ &
done
wait
for b in refGene knownGene wgEncodeGencodeCompVM2; do
  mkdir -p ../mouseBrowser_$release/liftoverbed/$b/C57B6J
  for f in browser_$release/liftoverbed/$b/C57B6J/*; do
    cp $f ../mouseBrowser_$release/liftoverbed/$b/C57B6J/ &
  done
done
wait
for b in refGene knownGene wgEncodeGencodeCompVM2; do
  cp browser_$release/liftoverbed/$b/C57B6J/C57B6J.as ../mouseBrowser_$release/liftoverbed/$b/C57B6J/ &
done
wait

echo '# Copying incidental files.'
# copy over files
for f in genomes.txt groups.txt hub.txt hubTreeModified.nw haltree.nw; do
  cp browser_$release/$f ../mouseBrowser_$release/ &
done
wait

echo '# Replacing 2bit symlinks with actual 2bit files.'
# replace 2bit symlinks with actual 2bit files
for g in $genomes; do
  rm -f ../mouseBrowser_$release/$g/$g.2bit;
  cp /hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_$release/$g.2bit ../mouseBrowser_$release/$g/$g.2bit &
done
wait

echo '# Fixing perimissions.'
# fix permissions
chmod -R g+w ../mouseBrowser_$release/* &> chmod_$release.log

echo '# Done.'
