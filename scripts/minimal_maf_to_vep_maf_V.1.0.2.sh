#!/bin/bash
# ======================================================================
# minimal_maf_to_vep_maf_V.1.0.2.sh
#
# 紐⑹쟻:
#   - minimal MAF (Chrom, Start, End, Ref, Alt, Sample ID留??덈뒗 ?곹깭)瑜?#      ?섑뵆蹂?VCF濡?履쇨컿 ??#      vcf2maf.pl 濡?annotation ?댁꽌
#      理쒖쥌?곸쑝濡?"?꾩꽦???쒖?) MAF" ?섎굹濡??⑹튂???뚯씠?꾨씪??
#
#   - 湲곕뒫:
#       * vcf2maf?먯꽌 annotation ?ㅽ뙣???섑뵆??蹂?대?
#         蹂꾨룄 ?뚯씪濡?紐⑥븘?????#         ??tmp_minimal_maf_to_vcf2maf/failed_variants_from_vcf2maf.tsv
#       * minimal MAF vs 理쒖쥌 vcf2maf MAF 瑜?鍮꾧탳?댁꽌
#         minimal?먮뒗 ?덉뿀?붾뜲 理쒖쥌 MAF?먮뒗 ?녿뒗 蹂??紐⑸줉 ?앹꽦
#         ??<OUTPUT_MAF>.dropped_from_minimal.tsv
#
# ?꾩젣:
#   - MSKCC vcf2maf ?ㅼ튂??#     ?? /path/mskcc-vcf2maf-754d68a/vcf2maf.pl
#   - hg19.fa (GRCh37) FASTA ?ㅼ튂??#   - Python3 + pandas ?ㅼ튂??# ======================================================================

set -euo pipefail

INPUT_MAF="${1:-}"
OUTPUT_MAF="${2:-}"

if [[ -z "${INPUT_MAF}" ]]; then
  echo "Usage: $0 <input_minimal_maf> [output_maf]" 1>&2
  exit 1
fi

if [[ ! -f "${INPUT_MAF}" ]]; then
  echo "[ERROR] Input MAF not found: ${INPUT_MAF}" 1>&2
  exit 1
fi

if [[ -z "${OUTPUT_MAF}" ]]; then
  base="${INPUT_MAF%.maf}"
  OUTPUT_MAF="${base}.vcf2maf.maf"
fi

VCF2MAF_DIR="/path/mskcc-vcf2maf-754d68a"
VCF2MAF_PL="${VCF2MAF_DIR}/vcf2maf.pl"

REF_FASTA="/path/.vep/hg19_genome/hg19.fa"
NCBI_BUILD="GRCh37"

CORES=16

TMP_ROOT="./tmp_minimal_maf_to_vcf2maf"
TMP_PY="${TMP_ROOT}/minimal_maf_split_to_vcf.py"
TMP_VCF_DIR="${TMP_ROOT}/vcfs"
TMP_MAF_DIR="${TMP_ROOT}/maf_from_vcf2maf"
VCF_LIST="${TMP_ROOT}/vcf_list.tsv"
FAILED_VARIANTS="${TMP_ROOT}/failed_variants_from_vcf2maf.tsv"

mkdir -p "${TMP_ROOT}" "${TMP_VCF_DIR}" "${TMP_MAF_DIR}"

echo -e "Sample_ID\tChromosome\tStart\tRef\tAlt\tReason" > "${FAILED_VARIANTS}"

cat > "${TMP_PY}" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd


def norm_chrom(chrom):
    """FASTA媛 chr1, chr2 ?뺤떇????VCF CHROM??留욎떠以?"""
    if chrom is None:
        return ""
    s = str(chrom).strip()
    if s == "":
        return ""
    if s.lower().startswith("chr"):
        return s
    return "chr" + s


def split_maf_to_vcfs(maf_in, out_dir, list_out):
    """
    minimal MAF瑜??섑뵆蹂?VCF濡?履쇨컿??

    ?꾩슂??MAF 而щ읆:
      - Tumor_Sample_Barcode
      - Chromosome
      - Start_Position
      - Reference_Allele
      - Tumor_Seq_Allele2
    VCF ?뺤떇:
      - #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample_id>
      - GT???쇰떒 0/1濡?紐⑤몢 ?ㅼ젙 (tumor-only, VAF ?뺣낫???놁쓬)
    """
    df = pd.read_csv(maf_in, sep="\t", dtype=str)
    required = [
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.stderr.write(f"[ERROR] MAF missing required columns: {', '.join(missing)}\n")
        sys.exit(1)

    df["_chr_norm"] = df["Chromosome"].apply(norm_chrom)
    df["_pos"] = df["Start_Position"]
    df["_ref"] = df["Reference_Allele"].fillna("").str.strip()
    df["_alt"] = df["Tumor_Seq_Allele2"].fillna("").str.strip()

    valid = (
        df["_chr_norm"].ne("")
        & df["_pos"].notna()
        & df["_ref"].ne("")
        & df["_alt"].ne("")
    )
    dfv = df.loc[valid].copy()

    if dfv.empty:
        sys.stderr.write("[ERROR] No valid rows in minimal MAF.\n")
        sys.exit(1)

    dfv["_pos"] = dfv["_pos"].astype(int)

    os.makedirs(out_dir, exist_ok=True)

    lines = []
    n_vcf = 0
    n_rows = 0

    for sid, sub in dfv.groupby("Tumor_Sample_Barcode"):
        vcf_path = os.path.join(out_dir, f"{sid}.from_minimal_maf.vcf")

        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=minimal_maf_to_vcf2maf\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sid)

            for _, r in sub.iterrows():
                chrom = r["_chr_norm"]
                pos = r["_pos"]
                ref = r["_ref"] if r["_ref"] != "" else "."
                alt = r["_alt"] if r["_alt"] != "" else "."

                hugo = r.get("Hugo_Symbol", "")
                info = f"ORIG_HUGO={hugo}" if isinstance(hugo, str) else "."

                line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1\n"
                f.write(line)
                n_rows += 1

        lines.append(f"{sid}\t{vcf_path}")
        n_vcf += 1

    with open(list_out, "w") as lf:
        for l in lines:
            lf.write(l + "\n")

    sys.stderr.write(
        f"[INFO] Wrote {n_vcf} per-sample VCFs with total {n_rows} variant rows.\n"
    )
    sys.stderr.write(f"[INFO] VCF list file: {list_out}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Split minimal MAF into per-sample VCFs for vcf2maf."
    )
    ap.add_argument("--maf-in", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--list-out", required=True)

    args = ap.parse_args()
    split_maf_to_vcfs(args.maf_in, args.out_dir, args.list_out)


if __name__ == "__main__":
    main()
PYEOF

chmod +x "${TMP_PY}"

echo "============================================"
echo "[1] minimal MAF ???섑뵆蹂?VCF ?앹꽦"
echo "============================================"

python3 "${TMP_PY}" \
  --maf-in "${INPUT_MAF}" \
  --out-dir "${TMP_VCF_DIR}" \
  --list-out "${VCF_LIST}"

echo ""
echo "============================================"
echo "[2] vcf2maf.pl ?ㅽ뻾 (?섑뵆蹂?VCF ???섑뵆蹂?MAF)"
echo "============================================"

if [[ ! -x "${VCF2MAF_PL}" ]]; then
  echo "[ERROR] vcf2maf.pl not found or not executable: ${VCF2MAF_PL}" 1>&2
  exit 1
fi

export FAILED_VARIANTS

if command -v parallel >/dev/null 2>&1; then
  echo "[INFO] GNU parallel 媛먯?????${CORES}肄붿뼱 蹂묐젹 ?ㅽ뻾"
  cat "${VCF_LIST}" | parallel -j "${CORES}" --colsep '\t' '
    sid={1}
    vcf={2}
    maf_out="'"${TMP_MAF_DIR}"'/${sid}.vcf2maf.maf"
    echo "   ??[${sid}] vcf2maf ?쒖옉"
    if perl "'"${VCF2MAF_PL}"'" \
      --input-vcf "$vcf" \
      --output-maf "$maf_out" \
      --ref-fasta "'"${REF_FASTA}"'" \
      --ncbi-build "'"${NCBI_BUILD}"'" \
      --tumor-id "$sid" \
      --vcf-tumor-id "$sid"; then
      echo "   ??[${sid}] ?꾨즺 ??$maf_out"
    else
      echo "   ??[${sid}] vcf2maf ?ㅽ뙣 ??${FAILED_VARIANTS}??湲곕줉"
      awk -v S="{1}" '\''BEGIN{FS="\t"} !/^#/ {print S"\t"$1"\t"$2"\t"$4"\t"$5"\tvcf2maf_failed"}'\'' "$vcf" >> "'"${FAILED_VARIANTS}"'"
    fi
  '
else
  echo "[INFO] GNU parallel ?놁쓬 ???쒖감 ?ㅽ뻾"
  while IFS=$'\t' read -r sid vcf; do
    maf_out="${TMP_MAF_DIR}/${sid}.vcf2maf.maf"
    echo "   ??[${sid}] vcf2maf ?쒖옉"
    if perl "${VCF2MAF_PL}" \
      --input-vcf "$vcf" \
      --output-maf "$maf_out" \
      --ref-fasta "${REF_FASTA}" \
      --ncbi-build "${NCBI_BUILD}" \
      --tumor-id "$sid" \
      --vcf-tumor-id "$sid"; then
      echo "   ??[${sid}] ?꾨즺 ??$maf_out"
    else
      echo "   ??[${sid}] vcf2maf ?ㅽ뙣 ??${FAILED_VARIANTS}??湲곕줉"
      awk -v S="$sid" 'BEGIN{FS="\t"} !/^#/ {print S"\t"$1"\t"$2"\t"$4"\t"$5"\tvcf2maf_failed"}' "$vcf" >> "${FAILED_VARIANTS}"
    fi
  done < "${VCF_LIST}"
fi

echo ""
echo "============================================"
echo "[3] ?섑뵆蹂?MAF瑜??섎굹??MAF濡??⑹튂湲?
echo "============================================"

shopt -s nullglob

maf_files=( "${TMP_MAF_DIR}"/*.vcf2maf.maf )

if [[ ${#maf_files[@]} -eq 0 ]]; then
  echo "[ERROR] No per-sample MAF files found in ${TMP_MAF_DIR}" 1>&2
  echo "[INFO] vcf2maf ?ㅽ뙣 蹂?대뒗 ${FAILED_VARIANTS} 瑜??뺤씤?섏꽭??"
  exit 1
fi

first_maf="${maf_files[0]}"

head -n 2 "${first_maf}" > "${OUTPUT_MAF}"
tail -n +3 "${first_maf}" >> "${OUTPUT_MAF}"

for f in "${maf_files[@]:1}"; do
  tail -n +3 "$f" >> "${OUTPUT_MAF}"
done

echo ""
echo "============================================"
echo "?꾨즺! vcf2maf annotated MAF ??${OUTPUT_MAF}"
echo " - Tumor_Sample_Barcode ???낅젰 minimal MAF??ID瑜?洹몃?濡??ъ슜"
echo " - 媛??섑뵆蹂?VCF 諛?以묎컙 MAF??${TMP_ROOT} ???덉쓬"
echo " - vcf2maf ?④퀎?먯꽌 ?ㅽ뙣??蹂?대뒗 ${FAILED_VARIANTS} ??湲곕줉??
echo "============================================"

echo ""
echo "============================================"
echo "[4] minimal MAF vs vcf2maf 理쒖쥌 MAF 鍮꾧탳 (drop??蹂??紐⑸줉 ?앹꽦)"
echo "============================================"

DIFF_PY="${TMP_ROOT}/compare_minimal_vs_vcf2maf.py"
DROPPED_OUT="${OUTPUT_MAF%.maf}.dropped_from_minimal.tsv"

cat > "${DIFF_PY}" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def norm_chr(x):
    s = str(x).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minimal", required=True)
    ap.add_argument("--final", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        df_min = pd.read_csv(args.minimal, sep="\t", dtype=str, comment="#")
    except Exception as e:
        sys.stderr.write(f"[ERROR] minimal MAF ?쎄린 ?ㅽ뙣: {e}\n")
        sys.exit(1)

    try:
        df_fin = pd.read_csv(args.final, sep="\t", dtype=str, comment="#")
    except Exception as e:
        sys.stderr.write(f"[ERROR] final MAF ?쎄린 ?ㅽ뙣: {e}\n")
        sys.exit(1)

    required = [
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    for c in required:
        if c not in df_min.columns:
            sys.stderr.write(f"[ERROR] minimal MAF??{c} 而щ읆???놁뒿?덈떎.\n")
            sys.exit(1)
        if c not in df_fin.columns:
            sys.stderr.write(f"[ERROR] final MAF??{c} 而щ읆???놁뒿?덈떎.\n")
            sys.exit(1)

    df_min["Chrom_norm"] = df_min["Chromosome"].apply(norm_chr)
    df_fin["Chrom_norm"] = df_fin["Chromosome"].apply(norm_chr)

    df_min["key"] = (
        df_min["Tumor_Sample_Barcode"].fillna("") + "|" +
        df_min["Chrom_norm"].fillna("") + "|" +
        df_min["Start_Position"].fillna("") + "|" +
        df_min["Reference_Allele"].fillna("") + "|" +
        df_min["Tumor_Seq_Allele2"].fillna("")
    )

    df_fin["key"] = (
        df_fin["Tumor_Sample_Barcode"].fillna("") + "|" +
        df_fin["Chrom_norm"].fillna("") + "|" +
        df_fin["Start_Position"].fillna("") + "|" +
        df_fin["Reference_Allele"].fillna("") + "|" +
        df_fin["Tumor_Seq_Allele2"].fillna("")
    )

    fin_keys = set(df_fin["key"])

    mask_drop = ~df_min["key"].isin(fin_keys)
    dropped = df_min.loc[mask_drop].copy()

    dropped.drop(columns=["Chrom_norm", "key"], inplace=True, errors="ignore")

    dropped.to_csv(args.out, sep="\t", index=False)

    sys.stderr.write(
        f"[DONE] minimal MAF 珥?{len(df_min)}媛?以?"
        f"理쒖쥌 MAF???녿뒗 蹂??{len(dropped)}媛쒕? {args.out} ?????n"
    )

if __name__ == "__main__":
    main()
PYEOF

chmod +x "${DIFF_PY}"

python3 "${DIFF_PY}" \
  --minimal "${INPUT_MAF}" \
  --final "${OUTPUT_MAF}" \
  --out "${DROPPED_OUT}"

echo ""
echo "============================================"
echo "minimal?먯꽌 drop??蹂??紐⑸줉: ${DROPPED_OUT}"
echo "============================================"
