# Fix snp_analysis/*/run_snp_analysis.sh
## Combine GATK, FreeBayes and Bcftools variant information: 
key_columns="#CHROM POS GENOTYPE GENE_ID GENE_NAME DESCRIPTION FEATURE_ID VARIANT_TYPE NU_CHANGE AA_CHANGE CDS_POS/CDS_LENGTH AA_POS/AA_LENGTH SAMPLE"
* Remove GENOTYPE from key_columns
Because different algorithmes may give different GENOTYPE values for the same variant
E.g.
#CHROM  POS     REF     ALT     QUAL    FORMAT  SAMPLE_INFO     ALLELE_FREQ     GENOTYPE        GENE_ID       GENE_NAME       DESCRIPTION     FEATURE_ID      VARIANT_TYPE    NU_CHANGE       AA_CHANGE    CDS_POS/CDS_LENGTH       AA_POS/AA_LENGTH        SAMPLE
## GATK
AE007317.1      161751  T       C       13396.1 GT:AD:DP:GQ:PL  1/1:2,300:302:99:13410,819,0    0.01,0.99     1/1     gene-spr0152-gene-spr0153       -|hk07  Conserved hypothetical protein|Histidine kinase       gene-spr0152-gene-spr0153       intergenic_region       n.161751T>C                          R6-WGT-9011-5
## FreeBayes
AE007317.1      161751  T       C       324.156 GT:DP:AD:RO:QR:AO:QA:GL 0/1:16:2,14:2:74:14:518:-42.146,0,-2.2094     0.12,0.88       0/1     gene-spr0152-gene-spr0153       -|hk07  Conserved hypothetical protein|Histidine kinase       gene-spr0152-gene-spr0153       intergenic_region       n.161751T>C  R6-WGT-9011-5

* Remove VARIANT_TYPE from key_columns
Because some variant type values are actually the same
E.g.
intergenic_region|synonymous_variant
synonymous_variant|intergenic_region
