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
## Differential Gene Expression
Gene ID outputs from htseq take several forms. Many genes include a "LOC" prefix as a placeholder until orthologs have been verified. LOC genes required conversion prior to database querying. More stable identifiers exist as lower case letters and digits such as "coq10b" and "pck1". These forms are ready for querying by most databases and do not require conversion.
After signficantly differentially expressed genes are identified (see [DESeq2 markdown](../)), they needed to be converted to a single format readable by a database, in this case the Maayan Lab [FishEnrichr](https://maayanlab.cloud/FishEnrichr/) database.
Conversion was done using the NCBI eFetch Utility and the following script
```bash
#!/bin/bash
STAMP="`date +'%FT%H%M'`"
FILENAME="multi-query_$STAMP.txt"
GENES="GENES_$FILENAME"
sp="/-\|"
sc=0
spin() {
   printf "\b${sp:sc++:1}"
   ((sc==${#sp})) && sc=0
}
endspin() {
   printf "\r%s\n" "$@"
}
cat next_batch.txt |
while read line
do
	spin
	if [[ ${line} =~ ^[a-zA-Z] ]]
	then
		echo "2." ${line} >> ./query_out/$FILENAME
	else
		efetch -db gene -id ${line} -format abstract >> ./query_out/$FILENAME
	fi
done
grep '^[0-9][.] ' ./query_out/$FILENAME >> ./query_out/$GENES
endspin
```
This script takes a set of genes with the "Loc-" or "gene-" prefix removed as an input in a file called next_batch.txt. LOC genes only contain numbers in this format and are fed into the eFetch gene db and the output in `-format abstract` is appended to an output file that is timestamped. Because this script would take 10 min. or more to run based on the number of input genes, there is also a status spinner so the user can see whether the query is still running after initiation. The output file includes unconverted genes preceded by "2." and converted genes with annotation preceded by "1."
```
1. prelid3b
PRELI domain containing 3 [Gadus morhua (Atlantic cod)]
Other Designations: PRELI domain containing protein 3B-like
Chromosome: 13
Annotation: Chromosome 13 NC_044060.1 (15022299..15027312)
ID: 115557485

2. atf3
```
The second output file created by the variable `$GENES` extracts just the gene names using the number and period at the beginning of the line to create a mostly clean gene list ready for database input.
```
2. sig_Tyr13_ordered_notrna_padj01
1. LOC115561805
2. b4galt7
1. LOC115559985
2. ppp1r10
```
In this example a name or identifier is included on the first line and maintained throughout. This list is then opened in excel using "." as a column delimiter to achieve clean lists for database input.

---
## Gene Ontology Annotation
There are several useful outputs from the FishEnrichr query. Annotation along different categories (this analysis focuses on "Biological Process") is available in a graphic and table format:

Bes17_finalish_padj01_GO_Biological_Process.png![image](https://user-images.githubusercontent.com/72388589/113826553-9669f080-9782-11eb-8544-5e6b80fb6486.png)

Term |	GO_id |	Overlap |	P-value	Adjusted | P-value |	Old P-value |	Old Adjusted P-value |	Z-score |	Combined Score |	Genes |
---|---|---|---|---|---|---|---|---|---|
inactivation of MAPK activity | GO:0000188 |	"4/6" |	7.78E-05 |	0.026178448 |	0.000192527 |	0.040969749 |	-4.038103465 |	38.20526364 |	dusp1;dusp2;dusp4;dusp5
hepatocyte differentiation |	GO:0070365 |	"3/6" |	0.002067377 |	0.183307423 |	0.00242818 |	0.123027786 |	-5.118114365 |	31.63749421 |	e2f8;apc;pck1
sister chromatid segregation |	GO:0000819 |	"6/16" | 6.94E-05 |	0.026178448 |	5.23E-05 |	0.020870389 |	-3.167958402 |	30.33427221 |	top2a;ncapd2;ncaph;mis12;ncapg;kif18a 

(table clipped for space, just an example)
From this table, significant GO terms can be extracted based on metric of choice and relevant biological pathway clusters can be constructed.

#### GO Visualization using simplifyEnrichment
*this analysis is on different data than the table above because that data was visualized manually*
```R
library(simplifyEnrichment)
library(magick)
IsFj_GO = read.delim(".../IsFj_padj01_finalish_GO_Biological_Process.txt", header = T)
go_sig = IsFj_GO[IsFj_GO$Adjusted.P.value<.05,]
go_id = go_sig$GO_id
mat = GO_similarity(go_id)
system.time({df = simplifyGO(mat,fontsize_range = c(10,18),draw_word_cloud = T)})
compare_clustering_methods(mat)
write.csv(df, ".../IsFj_GO_df.csv")
```

![IsFj_GO_cluster1](https://user-images.githubusercontent.com/72388589/113827614-cb2a7780-9783-11eb-8cda-400e34249125.png)
