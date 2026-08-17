# minimal_MAF_to_annotated_MAF_pipeline

Minimal-MAF conversion and re-annotation component of the CURE-NGS panel
harmonization framework.

> **Supported deployment:** use the unified
> [CURE-NGS Docker/OCI distribution](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework).
> This repository preserves the historical workflow. The container is
> published only from the umbrella repository, so this component may correctly
> display **No packages published**.

## Role in the unified project

| Item | Value |
| --- | --- |
| Historical responsibility | Split minimal MAF by sample, create VCF, annotate through vcf2maf |
| Supported commands | `cure-ngs minimal-maf-to-vcf` and `cure-ngs annotate-vcf` |
| Latest audited component release | `minimal_maf_to_vep_maf_V.1.0.2` |
| Default assembly | GRCh37/hg19 |

## Install the supported Docker distribution

1. Install [Docker Desktop](https://docs.docker.com/desktop/) or
   [Docker Engine](https://docs.docker.com/engine/install/).
2. Build both the small reviewer image and full annotation image:

```bash
git clone https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework.git
cd cure-ngs-panel-harmonization-framework
docker build --file docker/Dockerfile.core --tag cure-ngs-harmonizer:0.1.0-core .
docker build --file docker/Dockerfile --tag cure-ngs-harmonizer:0.1.0 .
```

After release `0.1.0` is published in the umbrella **Packages** panel:

```bash
docker pull ghcr.io/ncdcbioinformatics/cure-ngs-harmonizer:0.1.0-core
docker pull ghcr.io/ncdcbioinformatics/cure-ngs-harmonizer:0.1.0
```

## Verify and run this capability

Test minimal-MAF-to-VCF conversion without large reference downloads:

```bash
bash scripts/run_reviewer_demo.sh
```

Direct component command:

```bash
mkdir -p output
chmod 0777 output  # Linux: writable by the image's non-root UID 10001
docker run --rm \
  --volume "$PWD/examples:/examples:ro" \
  --volume "$PWD/output:/data/output" \
  cure-ngs-harmonizer:0.1.0-core minimal-maf-to-vcf \
  /examples/synthetic/minimal.grch37.maf /data/output/per-sample \
  --reference-fasta /examples/synthetic/tiny.grch37.fa \
  --assembly GRCh37
```

Full `annotate-vcf` execution uses the full image plus a real hg19 FASTA and
matching VEP 116 GRCh37 cache mounted read-only under `/references`.

## Historical standalone workflow

`minimal_maf_to_vep_maf_V.1.0.2.sh` remains available as release provenance.
It requires a manually installed vcf2maf checkout, Python, and hg19 FASTA. The
unified container removes those workstation-specific absolute-path assumptions.

## Documentation and test data

- [Project structure](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/PROJECT_STRUCTURE.md)
- [Minimal-MAF re-annotation commands](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/COMMAND_REFERENCE.md#minimal-maf-re-annotation-route)
- [Synthetic minimal MAF](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/examples/synthetic/minimal.grch37.maf)
- [Reference-data setup](https://github.com/NCDCbioinformatics/cure-ngs-panel-harmonization-framework/blob/main/docs/REFERENCE_DATA.md)

License: MIT. No CURE-NGS patient-level data are distributed here.
