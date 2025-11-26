# minimal_MAF_to_annotated_MAF_pipeline
Pipeline for creating annotated MAF from minimal MAF

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
   chmod +x minimal_maf_to_vep_maf_V.1.0.0.sh \
   ./minimal_maf_to_vep_maf_V.1.0.0.sh minimal_maf_from_hgvs_vep_V2.maf


# vcf2maf Execute permission required!!!!!
 chmod +x /path/mskcc-vcf2maf-754d68a/vcf2maf.pl
