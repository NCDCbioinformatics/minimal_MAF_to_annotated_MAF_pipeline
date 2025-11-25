# minimal_MAF_to_annotated_MAF_pipeline
minimal maf 에서 fully annotated MAF를 만드는 pipeline

# 목적:
   - minimal MAF (Chrom, Start, End, Ref, Alt, Sample ID만 있는 상태)를 \
     ▶ 샘플별 VCF로 쪼갠 뒤 \
     ▶ vcf2maf.pl 로 annotation 해서 \
     ▶ 최종적으로 "표준 MAF" 하나로 합치는 파이프라인

# 전제:
   - MSKCC vcf2maf 설치됨 \
     예) /home/wodn9614/mskcc-vcf2maf-754d68a/vcf2maf.pl \
   - hg19.fa (GRCh37) FASTA 설치됨 \
   - Python3 + pandas 설치됨

# 사용법:
   chmod +x minimal_maf_to_vcf2maf_maf.sh \
   ./minimal_maf_to_vcf2maf_maf.sh minimal_maf_from_hgvs_vep_V2.maf


#vcf2maf 실행 권한 필수!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 chmod +x /path/mskcc-vcf2maf-754d68a/vcf2maf.pl
