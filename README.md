# minimal_MAF_to_annotated_MAF_pipeline
Pipeline for creating annotated MAF from minimal MAF
<img width="2406" height="1335" alt="image" src="https://github.com/user-attachments/assets/40f2bca9-4c30-4d8c-b140-37978c3722a6" />

# purpose:
   - minimal MAF (only Chrom, Start, End, Ref, Alt, Sample ID) \
     ▶ After splitting into VCFs by sample \
     ▶ Annotate with vcf2maf.pl \
     ▶ Pipelines that ultimately merge into one "standard MAF"

# Requirements:
   - MSKCC vcf2maf installed \
     ex) /path/mskcc-vcf2maf-754d68a/vcf2maf.pl \
   - hg19.fa (GRCh37) FASTA installed \
   - Python3 + pandas installed

# How to use:
   chmod +x minimal_maf_to_vep_maf_V.1.0.2.sh \
   ./minimal_maf_to_vep_maf_V.1.0.2.sh minimal_maf_from_hgvs_vep_V2.maf


# vcf2maf Execute permission required!!!!!
 chmod +x /path/mskcc-vcf2maf-754d68a/vcf2maf.pl
