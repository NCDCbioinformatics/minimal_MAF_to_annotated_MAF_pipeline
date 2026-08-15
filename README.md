# minimal_MAF_to_annotated_MAF_pipeline

> **Reviewer and new-user deployment:** use the supported, version-pinned
> [CURE-NGS Docker/OCI distribution](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework#reviewer-quick-start).
> It replaces local absolute paths with bind-mounted input, output, and
> reference directories and validates the environment before annotation.

## Reproducible installation and test data

- [Clean-machine installation](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/INSTALLATION.md)
- [GRCh37 FASTA and release-matched VEP cache setup](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/REFERENCE_DATA.md)
- [Minimal-MAF re-annotation commands](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/COMMAND_REFERENCE.md#minimal-maf-re-annotation-route)
- [Network-free reviewer walkthrough](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/REVIEWER_REPRODUCTION.md)
- [Synthetic minimal MAF fixture](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/examples/synthetic/minimal.grch37.maf)

The latest audited historical release is
`minimal_maf_to_vep_maf_V.1.0.2`. Its tag, commit, asset size, and SHA-256 are
locked in the umbrella repository; the consolidated container is the supported
route for external reproduction.

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
## Publication context

This repository is a component of the CURE-NGS panel harmonization framework described in the manuscript "Multi-Institutional Harmonization Framework for Heterogeneous Panel-Based NGS in Precision Oncology."

Umbrella repository: https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework

## Software metadata

- Operating system(s): Linux; Windows users can run supported workflows via WSL where needed
- Programming language(s): Bash shell, Python
- License: MIT License
