## Pipeline from reads to gene ontology (GO) visualization
substantially borrowed from [Filipe Figueiredo](https://github.com/famfigueiredo/QuantSeq-January2020/blob/e8c15f31ae051790842a2ed34b26cc331d69428e/README.md)

### Demultiplexing with demuxFQ

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
# terminated because no files were being deposited in the demuxed files folder


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
# # Demux
#24.11.20
	# converted i5 to rc with added GT suffix based on i5 illumina guide and Alex's sample sheet

# lane 1 take 2
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take2.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take2.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_rc_gt_barcodes.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
	# terminated because no files were being deposited in the demuxed files folder


# lane 1 take 3 - new barcode file with only 3 columns, no space no extra identifiers
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take3.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take3.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_rc_gt_barcodes_new.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
	# summary says good, only 5% lost, but .fq.gz files are nowhere to be found. Only output in 01_lane1
	# is 4 .gz files corresponding to each of the i5s

# lane 1 take 4 - new barcode with 4 columns same as successful lane 2. sample sheet made in sublime and
# saved as .txt instead of excel export
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take4.lost.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane1-read1-take4.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane1_barcodes_sublime.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_1_FKDL202607271-1a_HF3N7CCX2_L7_1.fq.gz
	# started 16:13 25.11.20 finished before 21:00
	# BesFj1727_L1 indices came up as unexpected - need to trouble shoot, only 73 samples for L1 now
fastqc -o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/01_lane1/fastqc -t 60 -f fastq *.fq.gz



# lane 2 take 2
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/ \
-b /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane2-read1-take2.lostreads.fq.gz \
-s /data/ghi005/01_quantseq/02_unzipped-source-files/03_demux-summaries/lane2-read1-take2.summary.txt \
/data/ghi005/01_quantseq/03_useful-files/lane2_rc_gt_barcodes.txt \
/data/ghi005/01_quantseq/02_unzipped-source-files/01_fq.gz-files/GH_LANE_2_FKDL202607272-1a_HF3MCCCX2_L6_1.fq.gz 

# appears to be working, need to replace lane 1 master .fq.gz file because was overwritten to 0
```bash
fastqc -o /data/ghi005/01_quantseq/02_unzipped-source-files/02_demuxed-files/02_lane2/fastqc -t 60 -f fastq *.fq.gz
```
fastqc and multiqc done, fastqc files didn't end up in fastqc folder, need to play with that (how/where I run the actual script) (fixed)
