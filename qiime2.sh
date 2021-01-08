#!/bin/bash
#SBATCH -A p30777               # Allocation
#SBATCH -p long                # Queue
#SBATCH -t 10:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 4                    # Number of cores
#SBATCH --mem=0G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --mail-user=casar@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=ALL     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="qiime2"       # Name of job


#load singularity
module load singularity

#loop over directories 
for path in /projects/p30777/DeMMO_16s/DeMMO*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"


echo importing $dirname ...
##### import split, multiplexed paired-end reads from each run into .qza format #####
# input-path must be a folder containing the following files *exactly* named:
	# reverse.fastq.gz
	# forward.fastq.gz
	# barcodes.fastq.gz
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime tools import \
  --type EMPPairedEndSequences \
  --input-path ${path}/data \
  --output-path /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-emp-paired-end-sequences.qza'

echo demultiplexing $dirname ...
##### demultiplex sequence data ######
# provided metadata files (AKA mapping files) from Box folders must have the following features:
	# column "#SampleID" renamed to "sample-id"
	# column "BarcodeSequence" renamed to "barcode-sequence"
	# column "LinkerPrimerSequence" renamed to "linker-primer-sequence"
	# file must be saved with .tsv extension
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime demux emp-paired \
  --m-barcodes-file ${path}/$dirname'_map.txt' \
  --m-barcodes-column barcode-sequence \
  --i-seqs /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-emp-paired-end-sequences.qza' \
  --o-per-sample-sequences /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-demux.qza' \
  --o-error-correction-details /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-demux-details.qza' \
  --p-no-golay-error-correction # might need to include, per https://forum.qiime2.org/t/demux-with-very-little-reads/10161/11

echo picking ASVs for $dirname ...
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-demux.qza' \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-asv-table.qza' \
  --o-representative-sequences /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-rep-seqs.qza' \
  --o-denoising-stats /projects/p30777/Osburn2020/qiime2/$dirname/$dirname'-denoising-stats.qza'

  done

echo merging data ...
  ##### merge denoised data (only if you have data from multiple sequencing runs) #####
# merge table:
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime feature-table merge \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_Dec2015_Sep2016/DeMMO_Dec2015_Sep2016-asv-table.qza \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_May_to_Jul2016/DeMMO_May_to_Jul2016-asv-table.qza \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_Feb_to_Aug2017/DeMMO_Feb_to_Aug2017-asv-table.qza \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_Nov2017_to_April2018/DeMMO_Nov2017_to_April2018-asv-table.qza \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_Sep2018_Jun2019/DeMMO_Sep2018_Jun2019-asv-table.qza \
  --i-tables /projects/p30777/Osburn2020/qiime2/DeMMO_Dec2019/DeMMO_Dec2019-asv-table.qza \
  --o-merged-table /projects/p30777/Osburn2020/qiime2/Osburn2020-asv-table.qza

echo merging rep seqs ...
# merge representative-sequences:
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime feature-table merge-seqs \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_Dec2015_Sep2016/DeMMO_Dec2015_Sep2016-rep-seqs.qza \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_May_to_Jul2016/DeMMO_May_to_Jul2016-rep-seqs.qza \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_Feb_to_Aug2017/DeMMO_Feb_to_Aug2017-rep-seqs.qza \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_Nov2017_to_April2018/DeMMO_Nov2017_to_April2018-rep-seqs.qza \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_Sep2018_Jun2019/DeMMO_Sep2018_Jun2019-rep-seqs.qza \
  --i-data /projects/p30777/Osburn2020/qiime2/DeMMO_Dec2019/DeMMO_Dec2019-rep-seqs.qza \
  --o-merged-data /projects/p30777/Osburn2020/qiime2/Osburn2020-rep-seqs.qza

echo assigning taxonomy ...
 ##### assign taxonomy using pre-trained silva classifier for 515F/806R primer pair (via sklearn) #####
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime feature-classifier classify-sklearn \
  --i-classifier ~/silva-138-99-515-806-nb-classifier.qza \
  --i-reads /projects/p30777/Osburn2020/qiime2/Osburn2020-rep-seqs.qza \
  --p-reads-per-batch 10000 \
  --o-classification /projects/p30777/Osburn2020/qiime2/Osburn2020-taxonomy-Silva138.qza \
  --p-n-jobs 4

echo getting phylo tree ...
##### phylogenetic tree #####
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /projects/p30777/Osburn2020/qiime2/Osburn2020-rep-seqs.qza \
  --o-alignment /projects/p30777/Osburn2020/qiime2/Osburn2020-aligned-rep-seqs.qza \
  --o-masked-alignment /projects/p30777/Osburn2020/qiime2/Osburn2020-masked-aligned-rep-seqs.qza \
  --o-tree /projects/p30777/Osburn2020/qiime2/Osburn2020-unrooted-tree.qza \
  --o-rooted-tree /projects/p30777/Osburn2020/qiime2/Osburn2020-rooted-tree.qza

echo getting rarefaction curves ...
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime diversity alpha-rarefaction \
--i-table /projects/p30777/Osburn2020/qiime2/Osburn2020-asv-table.qza \
--i-phylogeny /projects/p30777/Osburn2020/qiime2/Osburn2020-rooted-tree.qza \
--p-max-depth 63000 \
--p-min-depth 500 \
--p-steps 500 \
--p-iterations 10 \
--m-metadata-file /projects/p30777/Osburn2020/qiime2/Osburn2020_metadata.txt \
--output-dir casar20201_rarefaction_data


echo rarefying data ...
##### rarefaction (sampling depth = 47468) #####
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime feature-table rarefy \
--i-table /projects/p30777/Osburn2020/qiime2/Osburn2020-asv-table.qza \
--p-sampling-depth 47468 \
--o-rarefied-table /projects/p30777/Osburn2020/qiime2/Casar2021-asv-table-rarefied-47468.qza


echo calculating diversity metrics ...

#calculate core metrics
singularity exec -B /projects/p30777 ~/qiime2_2020.8.sif qiime diversity core-metrics-phylogenetic \
--i-phylogeny /projects/p30777/Osburn2020/qiime2/Osburn2020-rooted-tree.qza \
--i-table /projects/p30777/Osburn2020/qiime2/Osburn2020-asv-table.qza \
--p-sampling-depth 47468 \
--m-metadata-file /projects/p30777/Osburn2020/qiime2/Osburn2020_metadata.txt \
--output-dir casar2021_core_metrics \
--p-n-jobs 16 \
--verbose \
&> core_metrics_samples.log

echo done!