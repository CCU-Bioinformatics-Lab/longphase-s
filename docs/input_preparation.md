## Input Preparation
- [Generate reference index](#generate-reference-index)
- [Generate alignment and index files](#generate-alignment-and-index-files)
- [Generate single nucleotide polymorphism (SNP) file](#generate-single-nucleotide-polymorphism-snp-file)
- [Generate Structural variation (SV)](#generate-structural-variation-sv-file)
- [Carry methylation tags to BAMs](#carry-methylation-tags-to-bams)

### Generate reference index
Index the reference genome with [samtools](https://github.com/samtools/samtools).
```
samtools faidx reference.fasta
```
### Generate alignment and index files
Produce read-to-reference alignment via [minimap2](https://github.com/lh3/minimap2) and sort/index the bam by [samtools](https://github.com/samtools/samtools).
```
# generate alignment flie with minimap2 according to the sequencing platform e.g. map-pb/map-ont/map-hifi
# Note that the MD-tag is required by sniffles (–MD).
minimap2 --MD -ax map-ont -t 10 reference.fasta reads.fastq -o alignment.sam

# sort alignment file
samtools sort -@ 10 alignment.sam -o alignment.bam

# index alignment file
samtools index -@ 10 alignment.bam
```
### Generate single nucleotide polymorphism (SNP) file
e.g. [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) or [Clair3](https://github.com/HKU-BAL/Clair3) pipeline.
```
INPUT_DIR={input data path}
OUTPUT_DIR={output data path}
BAM=alignment.bam
REF=reference.fasta
THREADS=10

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.7 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-t "${THREADS}" \
--ont_r9_guppy5_sup

# --ont_r9_guppy5_sup is preset for ONT R9.4.1 Guppy 5 "Sup" basecaller
# for ONT R10.4 Q20 reads: --ont_r10_q20
# for PacBio-HiFi reads: --hifi
```
### Generate Structural variation (SV) file
e.g. [sniffles](https://github.com/fritzsedlazeck/Sniffles) or [CuteSV](https://github.com/tjiangHIT/cuteSV).
```
# In sniffles1 please specofic --num_reads_report -1. For sniffles2 please specify --output-rnames instead.
sniffles -t 10 --num_reads_report -1 -m alignment.bam -v SV.vcf # for sniffles1
sniffles --threads 10 --output-rnames --input alignment.bam --vcf SV.vcf # for sniffles2

# cuteSV command for PacBio CLR data:
cuteSV alignment.bam reference.fasta SV.vcf work_dir --report_readid --genotype

# additional platform-specific parameters suggested by cuteSV
# PacBio CLR data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5
# PacBio CCS(HIFI) data: 
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
# ONT data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
```

### Carry methylation tags to BAMs
[The `-T` parameter](https://github.com/nanoporetech/dorado/issues/145) in samtools fastq extracts tags from the BAM file and stores them in the header of the FASTQ file. Please ensure that the BAM file includes both `MM` and `ML` tags and carried on in the following way.
```
samtools fastq -T '*' methylcall.raw.bam > methylcall.raw.fastq
```

Then, specify the `-y` option in minimap2 which appends tags stored in the FASTQ header into the BAM file.
```
minimap2 -ax map-ont -y reference.fasta methylcall.raw.fastq 
```