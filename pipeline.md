# Pipeline from reads to gene ontology (GO) visualization
substantially borrowed from [Filipe Figueiredo](https://github.com/famfigueiredo/QuantSeq-January2020/blob/e8c15f31ae051790842a2ed34b26cc331d69428e/README.md)

## Demultiplexing with demuxFQ

24.11.20
converted i5 to rc with added GT suffix based on i5 illumina guide and Alex's sample sheet

lane 1 take 2
```bash
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take2.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take2.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_rc_gt_barcodes.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
```
terminated because no files were being deposited in the demuxed files folder


lane 1 take 3 - new barcode file with only 3 columns, no space no extra identifiers
```bash
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take3.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take3.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_rc_gt_barcodes_new.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
```
summary says good, only 5% lost, but .fq.gz files are nowhere to be found. Only output in 01_lane1
is 4 .gz files corresponding to each of the i5s

lane 1 take 4 - new barcode with 4 columns same as successful lane 2. sample sheet made in sublime and
saved as .txt instead of excel export
```bash
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take4.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take4.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_barcodes_sublime.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
```
started 16:13 25.11.20 finished before 21:00
BesFj1727_L1 indices came up as unexpected - need to trouble shoot, only 73 samples for L1 now
```bash
fastqc -o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/fastqc -t 60 -f fastq *.fq.gz
```


lane 2 take 2
```bash
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane2-read1-take2.lostreads.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane2-read1-take2.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane2_rc_gt_barcodes.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_2_FKDL202607272-1a_HF3MCCCX2_L6_1.fq.gz 
```
appears to be working, need to replace lane 1 master .fq.gz file because was overwritten to 0

```bash
fastqc -o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/fastqc -t 60 -f fastq *.fq.gz
```
fastqc and multiqc done, fastqc files didn't end up in fastqc folder, need to play with that (how/where I run the actual script) (fixed)

---
## Adapter trimming with bbduk
Very sensitive script, ran into some copy paste issues that were eventually resolved creating the whole thing in nano.
Lane 1:
```bash
#!/bin/bash
cat /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/samplelist_nano_lane1.txt |
while read line
do
        bash /data/ghi005/bbmap/bbduk.sh  \
        in=/data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/$line.fq.gz \
        out=/data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/$line.clean.fq.gz \
        ref=/data/ghi005/01_quantseq/03_useful-files/polyA.fa,/data/ghi005/01_quantseq/03_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20

done
	# started 9:45 26/11/20
	# didn't work until I literally copied Filipe's working script to the lane 1 directory and just scratched out
	# the wrong lane numbers. Something really weird with my bbduk scripts in /03_useful-files/scripts

fastqc -o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/clean -t 60 -f fastq *.fq.gz
multiqc .
```
Lane 2:
```bash
#!/bin/bash
cat /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/samplelist_lane2.txt |
while read line
do
        bash /data/ghi005/bbmap/bbduk.sh  \
        in=/data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/$line.fq.gz \
        out=/data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/$line.clean.fq.gz \
        ref=/data/ghi005/01_quantseq/03_useful-files/polyA.fa,/data/ghi005/01_quantseq/03_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20

done
```
---
## Genome indexing with STAR
This step requires a reference genome including SIRV sequences (if used) in .gtf or .gff format
```bash
sudo STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /data/ghi005/01_quantseq/04_star_genomes/01_gadmorhua3.0-star \
--genomeSAindexNbases 13 \
--genomeFastaFiles /data/ghi005/01_quantseq/03_useful-files/genomes/gadMor3.0_SIRV-genomic.fasta \
--sjdbGTFfile /data/ghi005/01_quantseq/03_useful-files/genomes/gadMor3.0_SIRV-genomic.gtf |& tee genome-gen-out.txt
```
the output is a folder full of the components of the indexed genome, several GB in size
---
## Read alignment with STAR
This took several rounds, the the most useful output being achieved after completely relaxing the alignment parameters to be more accomodating of short reads.
The outcome of this "2short" optimized alignment was carried through the rest of analysis

```bash
#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/samplelist_nano_lane1.txt |
while read line
do
sudo STAR --runThreadN 30 --genomeDir /data/ghi005/01_quantseq/04_star-genomes/01_gadmorhua3.0-star --readFilesCommand zcat \
 --readFilesIn /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/clean/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/01_star-alignment-out/${line}_ |& tee starlog-lane1.txt \

done
```
average 82.8% aligned, finished 19:24
most unmapped are due to being "too short" - need to possibly relax this parameter to make sure we're not throwing away things
that are just a base or 2 too short

alignment 02_alignment2short
relaxed alignment parameters, everything associated with this second alignment will have "2short"
```bash
#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/samplelist_nano_lane1.txt |
while read line
do
sudo STAR --runThreadN 30 --genomeDir /data/ghi005/01_quantseq/04_star-genomes/01_gadmorhua3.0-star --readFilesCommand zcat \
 --readFilesIn /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/clean/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/01_star-alignment-out/02_alignment2short/${line}_ |& tee starlog-lane1.txt \

done
```
added --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
trying to match more of the "too short" reads
Nov 30 17:23:26 ..... started STAR run
Nov 30 18:53:02 ..... finished successfully

LANE 2

```bash
#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/samplelist_lane2.txt |
while read line
do
sudo STAR --runThreadN 30 --genomeDir /data/ghi005/01_quantseq/04_star-genomes/01_gadmorhua3.0-star --readFilesCommand zcat \
 --readFilesIn /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/clean/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/01_star-alignment-out/02_alignment2short/${line}_ |& tee starlog-lane2.txt \

done
```

Nov 28 09:09:58 ..... started STAR run
Nov 28 10:27:13 ..... finished successfully
added --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \
trying to match more of the "too short" reads
Dec 01 12:44:28 ..... started STAR run (short optimized)
---
