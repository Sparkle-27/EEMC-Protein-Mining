#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# EMCC protein filtering and construction of NR100, NR90 and NR50 sets
#
# Required input:
#   ROOT/proteins_batch1 ... ROOT/proteins_batch8
#   Each directory may contain *.faa or *.faa.gz files.
#
# Default filtering:
#   - length: 40-1000 aa, inclusive
#   - Prodigal completeness flag: partial=00
#   - valid protein alphabet
#
# Optional numeric completeness:
#   Export COMPLETENESS_TSV=/path/to/protein_completeness.tsv
#   Required columns: protein_id, completeness
#   The normalized protein_id format is GenomeID|OriginalProteinID.
#
# Independent clustering from the same filtered input:
#   NR100: identity 1.00, coverage 1.00
#   NR90 : identity 0.90, coverage 0.90
#   NR50 : identity 0.50, coverage 0.90
#
# Usage:
#   bash 02_run_EMCC_filter_cluster.sh /path/to/EMCC_root
#
# Background usage:
#   nohup env THREADS=32 bash 02_run_EMCC_filter_cluster.sh \
#     /path/to/EMCC_root > EMCC_filter_cluster.nohup.log 2>&1 &
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOFTWARE_CONFIG="${SOFTWARE_CONFIG:-${SCRIPT_DIR}/00_software_paths.sh}"
FILTER_PROGRAM="${FILTER_PROGRAM:-${SCRIPT_DIR}/01_filter_EMCC_proteins.py}"

[[ -r "$SOFTWARE_CONFIG" ]] || {
    echo "ERROR: software configuration not found: $SOFTWARE_CONFIG" >&2
    exit 1
}
# shellcheck source=/dev/null
source "$SOFTWARE_CONFIG"
check_emcc_software
print_emcc_software

ROOT="${1:-$PWD}"
OUTDIR="${2:-${ROOT}/EMCC_protein_filtering}"

THREADS="${THREADS:-32}"
MIN_LENGTH="${MIN_LENGTH:-40}"
MAX_LENGTH="${MAX_LENGTH:-1000}"
REQUIRE_PARTIAL00="${REQUIRE_PARTIAL00:-1}"

MIN_COMPLETENESS="${MIN_COMPLETENESS:-90}"
COMPLETENESS_TSV="${COMPLETENESS_TSV:-}"

IDENTITY_NR100="${IDENTITY_NR100:-1.00}"
COVERAGE_NR100="${COVERAGE_NR100:-1.00}"

IDENTITY_NR90="${IDENTITY_NR90:-0.90}"
COVERAGE_NR90="${COVERAGE_NR90:-0.90}"

IDENTITY_NR50="${IDENTITY_NR50:-0.50}"
COVERAGE_NR50="${COVERAGE_NR50:-0.90}"

NR100_WORKFLOW="${NR100_WORKFLOW:-easy-linclust}"
NR90_WORKFLOW="${NR90_WORKFLOW:-easy-linclust}"
NR50_WORKFLOW="${NR50_WORKFLOW:-easy-cluster}"

TMP_BASE="${TMP_BASE:-${OUTDIR}/tmp}"
KEEP_MMSEQS_ALL_SEQS="${KEEP_MMSEQS_ALL_SEQS:-0}"

mkdir -p \
    "${OUTDIR}/00_metadata" \
    "${OUTDIR}/01_filtered" \
    "${OUTDIR}/02_NR100" \
    "${OUTDIR}/03_NR90" \
    "${OUTDIR}/04_NR50" \
    "${OUTDIR}/05_summary" \
    "${OUTDIR}/logs" \
    "$TMP_BASE"

LOG="${OUTDIR}/logs/workflow.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

validate_fraction() {
    local name="$1"
    local value="$2"
    "$PYTHON3" - "$name" "$value" <<'PY'
import sys
name, raw = sys.argv[1], sys.argv[2]
try:
    value = float(raw)
except ValueError:
    raise SystemExit(f"{name} is not numeric: {raw}")
if not (0.0 < value <= 1.0):
    raise SystemExit(f"{name} must be in (0,1]: {value}")
PY
}

for value in "$MIN_LENGTH" "$MAX_LENGTH" "$THREADS"; do
    [[ "$value" =~ ^[0-9]+$ ]] || die "Length/thread parameters must be integers."
done
(( MIN_LENGTH >= 1 && MAX_LENGTH >= MIN_LENGTH )) || die "Invalid length range."

validate_fraction IDENTITY_NR100 "$IDENTITY_NR100"
validate_fraction COVERAGE_NR100 "$COVERAGE_NR100"
validate_fraction IDENTITY_NR90 "$IDENTITY_NR90"
validate_fraction COVERAGE_NR90 "$COVERAGE_NR90"
validate_fraction IDENTITY_NR50 "$IDENTITY_NR50"
validate_fraction COVERAGE_NR50 "$COVERAGE_NR50"

for workflow in "$NR100_WORKFLOW" "$NR90_WORKFLOW" "$NR50_WORKFLOW"; do
    case "$workflow" in
        easy-cluster|easy-linclust) ;;
        *) die "Unsupported MMseqs2 workflow: $workflow" ;;
    esac
done

[[ -r "$FILTER_PROGRAM" ]] || die "Filter program not found: $FILTER_PROGRAM"
[[ -d "$ROOT" ]] || die "Input root directory not found: $ROOT"

log "Input root: $ROOT"
log "Output directory: $OUTDIR"
log "Length filter: ${MIN_LENGTH}-${MAX_LENGTH} aa"
log "Require Prodigal partial=00: $REQUIRE_PARTIAL00"
if [[ -n "$COMPLETENESS_TSV" ]]; then
    [[ -r "$COMPLETENESS_TSV" ]] || die "Completeness TSV not readable: $COMPLETENESS_TSV"
    log "External numeric completeness filter: >=${MIN_COMPLETENESS}%"
    log "Completeness table: $COMPLETENESS_TSV"
else
    log "External numeric completeness table: not supplied"
    log "Completeness is operationalized as Prodigal partial=00 only."
fi

# ---------- Record parameters and software ----------
write_emcc_software_versions "${OUTDIR}/00_metadata/software_versions.tsv"

{
    printf "parameter\tvalue\tdescription\n"
    printf "ROOT\t%s\tinput project root\n" "$ROOT"
    printf "OUTDIR\t%s\tworkflow output directory\n" "$OUTDIR"
    printf "MIN_LENGTH\t%s\tinclusive minimum protein length\n" "$MIN_LENGTH"
    printf "MAX_LENGTH\t%s\tinclusive maximum protein length\n" "$MAX_LENGTH"
    printf "REQUIRE_PARTIAL00\t%s\tkeep only Prodigal partial=00\n" "$REQUIRE_PARTIAL00"
    printf "COMPLETENESS_TSV\t%s\toptional external numeric completeness table\n" \
        "${COMPLETENESS_TSV:-NOT_SUPPLIED}"
    printf "MIN_COMPLETENESS\t%s\texternal completeness threshold if table supplied\n" \
        "$MIN_COMPLETENESS"
    printf "IDENTITY_NR100\t%s\tMMseqs2 minimum sequence identity\n" "$IDENTITY_NR100"
    printf "COVERAGE_NR100\t%s\tMMseqs2 coverage of longer sequence\n" "$COVERAGE_NR100"
    printf "IDENTITY_NR90\t%s\tMMseqs2 minimum sequence identity\n" "$IDENTITY_NR90"
    printf "COVERAGE_NR90\t%s\tMMseqs2 coverage of longer sequence\n" "$COVERAGE_NR90"
    printf "IDENTITY_NR50\t%s\tMMseqs2 minimum sequence identity\n" "$IDENTITY_NR50"
    printf "COVERAGE_NR50\t%s\tMMseqs2 coverage of longer sequence\n" "$COVERAGE_NR50"
    printf "COV_MODE\t0\taligned residues divided by longer sequence length\n"
    printf "ALIGNMENT_MODE\t3\tidentity from identical aligned residues/alignment columns\n"
    printf "THREADS\t%s\tMMseqs2 threads\n" "$THREADS"
    printf "NR100_WORKFLOW\t%s\tMMseqs2 easy workflow\n" "$NR100_WORKFLOW"
    printf "NR90_WORKFLOW\t%s\tMMseqs2 easy workflow\n" "$NR90_WORKFLOW"
    printf "NR50_WORKFLOW\t%s\tMMseqs2 easy workflow\n" "$NR50_WORKFLOW"
} > "${OUTDIR}/00_metadata/run_parameters.tsv"

# ---------- Discover input files ----------
FILE_LIST="${OUTDIR}/00_metadata/input_faa_files.txt"
: > "$FILE_LIST"

for batch in {1..8}; do
    batch_dir="${ROOT}/proteins_batch${batch}"
    [[ -d "$batch_dir" ]] || die "Expected batch directory is missing: $batch_dir"

    find "$batch_dir" -type f \( \
        -name '*.faa' -o -name '*.faa.gz' -o \
        -name '*.fasta' -o -name '*.fasta.gz' -o \
        -name '*.fa' -o -name '*.fa.gz' \
    \) -print >> "$FILE_LIST"
done

sort -u "$FILE_LIST" -o "$FILE_LIST"
INPUT_FILE_COUNT="$(wc -l < "$FILE_LIST" | tr -d ' ')"
(( INPUT_FILE_COUNT > 0 )) || die "No protein FASTA files found in proteins_batch1-8."
log "Protein FASTA files detected: $INPUT_FILE_COUNT"

# FAA basenames are used as genome IDs and must therefore be unique.
DUPLICATE_GENOMES="${OUTDIR}/00_metadata/duplicate_genome_ids.txt"
"$PYTHON3" - "$FILE_LIST" "$DUPLICATE_GENOMES" <<'PY'
import sys
from collections import Counter
from pathlib import Path

file_list = Path(sys.argv[1])
out = Path(sys.argv[2])

def genome_id(path):
    name = Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in (".faa", ".fasta", ".fa"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return name

ids = [genome_id(x) for x in file_list.read_text().splitlines() if x.strip()]
duplicates = sorted(k for k, v in Counter(ids).items() if v > 1)
out.write_text("".join(x + "\n" for x in duplicates))
PY

if [[ -s "$DUPLICATE_GENOMES" ]]; then
    cat "$DUPLICATE_GENOMES" >&2
    die "Duplicate genome IDs detected. Rename repeated FASTA basenames."
fi

# ---------- Filter and normalize ----------
FILTERED_FASTA="${OUTDIR}/01_filtered/EMCC.filtered.L${MIN_LENGTH}-${MAX_LENGTH}.partial00.faa"
FILTERED_MANIFEST="${OUTDIR}/01_filtered/retained_proteins.tsv"
DISCARDED_MANIFEST="${OUTDIR}/01_filtered/discarded_proteins.tsv"
FILTER_SUMMARY="${OUTDIR}/05_summary/filter_summary.tsv"

filter_command=(
    "$PYTHON3" "$FILTER_PROGRAM"
    --file-list "$FILE_LIST"
    --output-fasta "$FILTERED_FASTA"
    --kept-table "$FILTERED_MANIFEST"
    --discarded-table "$DISCARDED_MANIFEST"
    --summary "$FILTER_SUMMARY"
    --min-length "$MIN_LENGTH"
    --max-length "$MAX_LENGTH"
)

if [[ "$REQUIRE_PARTIAL00" == "1" ]]; then
    filter_command+=(--require-partial00)
fi

if [[ -n "$COMPLETENESS_TSV" ]]; then
    filter_command+=(
        --completeness-table "$COMPLETENESS_TSV"
        --min-completeness "$MIN_COMPLETENESS"
    )
fi

log "Filtering and normalizing protein sequences."
printf '%q ' "${filter_command[@]}" > "${OUTDIR}/00_metadata/filter_command.txt"
printf '\n' >> "${OUTDIR}/00_metadata/filter_command.txt"
"${filter_command[@]}"

[[ -s "$FILTERED_FASTA" ]] || die "Filtered FASTA is empty."

# ---------- MMseqs2 clustering ----------
run_mmseqs_cluster() {
    local label="$1"
    local identity="$2"
    local coverage="$3"
    local workflow="$4"
    local result_dir="$5"

    local prefix="${result_dir}/EMCC_${label}"
    local tmp_dir="${TMP_BASE}/${label}"

    mkdir -p "$result_dir"
    rm -rf "$tmp_dir"
    mkdir -p "$tmp_dir"

    log "Building ${label}: identity=${identity}; coverage=${coverage}; workflow=${workflow}"

    local cmd=(
        "$MMSEQS" "$workflow"
        "$FILTERED_FASTA"
        "$prefix"
        "$tmp_dir"
        --min-seq-id "$identity"
        -c "$coverage"
        --cov-mode 0
        --alignment-mode 3
        --threads "$THREADS"
        --remove-tmp-files 1
    )

    printf '%q ' "${cmd[@]}" > "${OUTDIR}/00_metadata/${label}_command.txt"
    printf '\n' >> "${OUTDIR}/00_metadata/${label}_command.txt"

    "${cmd[@]}"

    [[ -s "${prefix}_rep_seq.fasta" ]] || die "${label} representative FASTA missing."
    [[ -s "${prefix}_cluster.tsv" ]] || die "${label} cluster mapping missing."

    cp "${prefix}_rep_seq.fasta" "${result_dir}/${label}.representatives.faa"
    cp "${prefix}_cluster.tsv" "${result_dir}/${label}.representative_to_member.tsv"

    if [[ "$KEEP_MMSEQS_ALL_SEQS" != "1" ]]; then
        rm -f "${prefix}_all_seqs.fasta"
    fi
}

# The three sets are independently generated from the same filtered FASTA.
run_mmseqs_cluster \
    NR100 "$IDENTITY_NR100" "$COVERAGE_NR100" "$NR100_WORKFLOW" \
    "${OUTDIR}/02_NR100"

run_mmseqs_cluster \
    NR90 "$IDENTITY_NR90" "$COVERAGE_NR90" "$NR90_WORKFLOW" \
    "${OUTDIR}/03_NR90"

run_mmseqs_cluster \
    NR50 "$IDENTITY_NR50" "$COVERAGE_NR50" "$NR50_WORKFLOW" \
    "${OUTDIR}/04_NR50"

# ---------- Count and summarize ----------
count_fasta() {
    "$PYTHON3" - "$1" <<'PY'
import gzip
import sys
from pathlib import Path

path = Path(sys.argv[1])
opener = gzip.open if path.name.endswith(".gz") else open
n = 0
with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith(">"):
            n += 1
print(n)
PY
}

N_FILTERED="$(count_fasta "$FILTERED_FASTA")"
N_NR100="$(count_fasta "${OUTDIR}/02_NR100/NR100.representatives.faa")"
N_NR90="$(count_fasta "${OUTDIR}/03_NR90/NR90.representatives.faa")"
N_NR50="$(count_fasta "${OUTDIR}/04_NR50/NR50.representatives.faa")"

CLUSTER_SUMMARY="${OUTDIR}/05_summary/cluster_set_summary.tsv"
"$PYTHON3" - \
    "$N_FILTERED" "$N_NR100" "$N_NR90" "$N_NR50" \
    "$IDENTITY_NR100" "$COVERAGE_NR100" \
    "$IDENTITY_NR90" "$COVERAGE_NR90" \
    "$IDENTITY_NR50" "$COVERAGE_NR50" \
    > "$CLUSTER_SUMMARY" <<'PY'
import sys

(
    n_filtered, n100, n90, n50,
    id100, cov100, id90, cov90, id50, cov50
) = sys.argv[1:]

n_filtered = int(n_filtered)
rows = [
    ("FILTERED", "NA", "NA", n_filtered),
    ("NR100", id100, cov100, int(n100)),
    ("NR90", id90, cov90, int(n90)),
    ("NR50", id50, cov50, int(n50)),
]

print("set\tidentity_cutoff\tcoverage_cutoff\tcoverage_mode\tprotein_count\tretention_vs_filtered")
for name, identity, coverage, count in rows:
    retention = count / n_filtered if n_filtered else 0.0
    cov_mode = "NA" if name == "FILTERED" else "0_longer_sequence"
    print(
        f"{name}\t{identity}\t{coverage}\t{cov_mode}\t"
        f"{count}\t{retention:.6f}"
    )
PY

log "Filtered proteins: $N_FILTERED"
log "NR100 representatives: $N_NR100"
log "NR90 representatives:  $N_NR90"
log "NR50 representatives:  $N_NR50"

if (( N_NR100 > N_FILTERED || N_NR90 > N_FILTERED || N_NR50 > N_FILTERED )); then
    die "A clustered set is larger than the filtered input; inspect the run."
fi

if ! (( N_NR100 >= N_NR90 && N_NR90 >= N_NR50 )); then
    log "WARNING: counts are not strictly monotonic. Inspect commands, version and logs."
fi

cat > "${OUTDIR}/05_summary/README.txt" <<EOF
EMCC protein filtering and three-level clustering completed.

PURPOSE
  Filter Prodigal proteins and construct three independently clustered
  representative protein sets for downstream annotation and structure mining.

INPUT
  ${ROOT}/proteins_batch1 ... proteins_batch8
  Temporary directories proteins_tmp1 and proteins_tmp2 are not used.

FILTERING
  Protein length: ${MIN_LENGTH}-${MAX_LENGTH} aa, inclusive
  Prodigal full-CDS flag required: partial=00 = ${REQUIRE_PARTIAL00}
  External numeric completeness table:
    ${COMPLETENESS_TSV:-not supplied}
  External completeness threshold:
    ${MIN_COMPLETENESS}% (only applied when a table is supplied)

CLUSTERING
  NR100: identity ${IDENTITY_NR100}; coverage ${COVERAGE_NR100}
  NR90 : identity ${IDENTITY_NR90}; coverage ${COVERAGE_NR90}
  NR50 : identity ${IDENTITY_NR50}; coverage ${COVERAGE_NR50}
  MMseqs2 coverage mode: 0, aligned residues / longer sequence length
  MMseqs2 alignment mode: 3, exact aligned-residue identity calculation

KEY OUTPUTS
  Filtered full set:
    ${FILTERED_FASTA}

  NR100 representatives:
    ${OUTDIR}/02_NR100/NR100.representatives.faa
  NR100 mapping:
    ${OUTDIR}/02_NR100/NR100.representative_to_member.tsv

  NR90 representatives:
    ${OUTDIR}/03_NR90/NR90.representatives.faa
  NR90 mapping:
    ${OUTDIR}/03_NR90/NR90.representative_to_member.tsv

  NR50 representatives:
    ${OUTDIR}/04_NR50/NR50.representatives.faa
  NR50 mapping:
    ${OUTDIR}/04_NR50/NR50.representative_to_member.tsv

  Filter summary:
    ${FILTER_SUMMARY}
  Cluster summary:
    ${CLUSTER_SUMMARY}

INTERPRETATION
  NR100 removes exact full-length duplicates.
  NR90 is recommended for routine functional annotation and candidate screening.
  NR50 is recommended for broader protein-family and structural-space analyses.
  Cluster annotations must be mapped back to the original proteins using the
  representative_to_member tables.
EOF

log "Workflow completed successfully."
log "Run log: $LOG"
log "Summary: $CLUSTER_SUMMARY"
