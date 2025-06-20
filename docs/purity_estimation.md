## Tumor purity estimation command
This command using tumor-normal pair BAM and VCF files, along with haplotype information, and outputs the prediction file.
```
longphase-s estimate_purity \
-s phased_normal_snp.vcf \
-b normal.bam \
--tumor-snp-file tumor_snp.vcf \
--tumor-bam-file tumor.bam \
-r reference.fasta \
-t 8 \
-o output_prefix
```

### The complete list of purity estimation parameters

```
Usage:  estimate_purity [OPTION] ... READSFILE
      --help                          display this help and exit.

required arguments:
      -s, --snp-file=NAME             input phased normal sample SNP VCF file.
      -b, --bam-file=NAME             input normal sample BAM file.
      --tumor-snp-file=NAME           input tumor sample SNP VCF file.
      --tumor-bam-file=NAME           input tumor sample BAM file.
      -r, --reference=NAME            reference FASTA.

optional arguments:
      --tagSupplementary              include supplementary alignments in haplotype statistics. default:true
      -q, --qualityThreshold=Num      exclude reads with mapping quality below threshold. default:20
      -p, --percentageThreshold=Num   include alignments in statistics based on haplotype allele percentage.
                                      alignments without clear haplotype assignment are excluded. default:0.6
      -t, --threads=Num               number of thread. default:1
      -o, --out-prefix=NAME           prefix of tumor purity estimation result. default:result
      --region=REGION                 tumor purity estimation include only reads/variants overlapping those regions. default:""(all regions)
                                      input format:chrom (consider entire chromosome)
                                                   chrom:start (consider region from this start to end of chromosome)
                                                   chrom:start-end
```