#!/usr/bin/env bash
# End-to-end T2tt jet-bin-fix comparison pipeline:
#   .result files -> datacards -> combine -> comparison plot + 174-bin plot
#
# Run this from the HOST shell (NOT already inside the OSCA/coffea container),
# from the repository root. It enters/exits the container itself as needed:
#   scripts/combine/run_full_pipeline.sh
#
# Defaults point at the full-statistics condor rerun output directories.
# Pass --test to use the quick local 10k-event result files instead.
#
# Optional overrides (env vars):
#   SIGNAL_DIR      default: analysis_products/old_results/26-09-10_hadsusy_signal_full_allmass
#                            (--test: analysis_products/old_results/26-09-10_hadsusy_signal_10k)
#   BACKGROUND_DIR  default: analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed
#                            (--test: analysis_products/old_results/26-09-10_hadsusy_backgrounds_10k)
#   COMBINE_DIR     default: analysis_products/combine/split_signals_jetbins_fixed
#   SEARCHBIN_PLOT_DIR   default: analysis_products/plots/diagnostics/search_bin_diagnostics
#   SEARCHBIN_COMBINED_OUT default: analysis_products/plots/diagnostics/search_bin_diagnostics_combined.png
#                          (primary plot: all mLSP signals summed into one "signal" line)
#   UV_PROJECT_ENVIRONMENT  default: $PWD/.venv-local10k

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

TEST_MODE=false
for arg in "$@"; do
    if [[ "$arg" == "--test" ]]; then
        TEST_MODE=true
    fi
done

if [[ "$TEST_MODE" == "true" ]]; then
    DEFAULT_SIGNAL_DIR="analysis_products/old_results/26-09-10_hadsusy_signal_10k"
    DEFAULT_BACKGROUND_DIR="analysis_products/old_results/26-09-10_hadsusy_backgrounds_10k"
else
    DEFAULT_SIGNAL_DIR="analysis_products/old_results/26-09-10_hadsusy_signal_full_allmass"
    DEFAULT_BACKGROUND_DIR="analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed"
fi

SIGNAL_DIR="${SIGNAL_DIR:-$DEFAULT_SIGNAL_DIR}"
BACKGROUND_DIR="${BACKGROUND_DIR:-$DEFAULT_BACKGROUND_DIR}"
MERGED_DIR="${MERGED_DIR:-.application_data/pipeline_results_merged}"
COMBINE_DIR="${COMBINE_DIR:-analysis_products/combine/split_signals_jetbins_fixed}"
LIMITS_CSV="${LIMITS_CSV:-$COMBINE_DIR/limits.csv}"
COMPARISON_OUT="${COMPARISON_OUT:-$COMBINE_DIR/t2tt_hepdata_comparison}"
SEARCHBIN_PLOT_DIR="${SEARCHBIN_PLOT_DIR:-analysis_products/plots/diagnostics/search_bin_diagnostics}"
SEARCHBIN_COMBINED_OUT="${SEARCHBIN_COMBINED_OUT:-analysis_products/plots/diagnostics/search_bin_diagnostics_combined.png}"
UV_PROJECT_ENVIRONMENT="${UV_PROJECT_ENVIRONMENT:-$PWD/.venv-local10k}"
ALLOW_MISSING_FLAG=""
if [[ "$TEST_MODE" == "true" ]]; then
    ALLOW_MISSING_FLAG="--allow-missing-backgrounds"
fi

echo "Mode: $([[ "$TEST_MODE" == "true" ]] && echo '--test (10k local)' || echo 'full condor rerun')"
echo "SIGNAL_DIR=$SIGNAL_DIR"
echo "BACKGROUND_DIR=$BACKGROUND_DIR"

echo "=== 0/4: merging signal + background result files into $MERGED_DIR ==="
rm -rf "$MERGED_DIR"
mkdir -p "$MERGED_DIR"
ln -sf "$(realpath "$SIGNAL_DIR")"/*.result "$MERGED_DIR"/
ln -sf "$(realpath "$BACKGROUND_DIR")"/*.result "$MERGED_DIR"/

echo "=== 1/4: generating + running per-signal datacards (inside container) ==="
UV_PROJECT_ENVIRONMENT="$UV_PROJECT_ENVIRONMENT" MERGED_DIR="$MERGED_DIR" COMBINE_DIR="$COMBINE_DIR" ALLOW_MISSING_FLAG="$ALLOW_MISSING_FLAG" ./setup.sh <<'EOF'
set -euo pipefail
uv run --no-sync python3 scripts/combine/make_split_signal_datacards.py \
    "$MERGED_DIR" --output-dir "$COMBINE_DIR" --run --overwrite $ALLOW_MISSING_FLAG
EOF

echo "=== 2/4: running combine (host apptainer) ==="
python3 scripts/combine/run_split_signal_combine.py "$COMBINE_DIR" --keep-going

echo "=== 3/4: collecting limits + plotting comparison (inside container) ==="
UV_PROJECT_ENVIRONMENT="$UV_PROJECT_ENVIRONMENT" COMBINE_DIR="$COMBINE_DIR" LIMITS_CSV="$LIMITS_CSV" COMPARISON_OUT="$COMPARISON_OUT" ./setup.sh <<'EOF'
set -euo pipefail
uv run --no-sync python3 scripts/combine/collect_split_signal_limits.py \
    "$COMBINE_DIR" --output "$LIMITS_CSV"
uv run --no-sync python3 scripts/combine/plot_t2tt_hepdata_comparison.py \
    --limits-csv "$LIMITS_CSV" --output "$COMPARISON_OUT"
EOF

echo "=== 4/4: 174-bin search-region diagnostic plots (inside container) ==="
UV_PROJECT_ENVIRONMENT="$UV_PROJECT_ENVIRONMENT" COMBINE_DIR="$COMBINE_DIR" SEARCHBIN_PLOT_DIR="$SEARCHBIN_PLOT_DIR" SEARCHBIN_COMBINED_OUT="$SEARCHBIN_COMBINED_OUT" ./setup.sh <<'EOF'
set -euo pipefail
mkdir -p "$SEARCHBIN_PLOT_DIR" "$(dirname "$SEARCHBIN_COMBINED_OUT")"

plot_diag() {
    local shapes_file="$1" output="$2" signal_hist="$3"
    local background_args
    background_args=$(uv run --no-sync python3 -c "
import uproot
keys = {k.split(';')[0] for k in uproot.open('$shapes_file').keys()}
defaults = ['qcd_inclusive_2024','tt_hadronic_2024','tt_semileptonic_2024','wjets_2024','zjets_2024']
print(' '.join(f'--background {b}' for b in defaults if b in keys))
")
    uv run --no-sync python3 scripts/combine/plot_search_bin_diagnostics.py \
        "$shapes_file" --signal "$signal_hist" $background_args -o "$output"
}

# One plot per mass point.
shapes_files=()
for point_dir in "$COMBINE_DIR"/T2tt_mStop-*_mLSP-*; do
    point=$(basename "$point_dir")
    shapes_file="$point_dir/shapes_bins_hadsusy.root"
    [[ -f "$shapes_file" ]] || continue
    shapes_files+=("$shapes_file")
    signal_hist=$(uv run --no-sync python3 -c "
import uproot
keys = {k.split(';')[0] for k in uproot.open('$shapes_file').keys()}
print(next(k for k in keys if k.startswith('2stop_')))
")
    plot_diag "$shapes_file" "$SEARCHBIN_PLOT_DIR/${point}.png" "$signal_hist"
done

# Primary combined plot: all mLSP signals summed into one "signal" line.
COMBINED_SHAPES="$COMBINE_DIR/shapes_all_signals_combined.root"
uv run --no-sync python3 scripts/combine/build_combined_signal_shapes.py \
    "${shapes_files[@]}" --signal-histogram "$signal_hist" -o "$COMBINED_SHAPES"
plot_diag "$COMBINED_SHAPES" "$SEARCHBIN_COMBINED_OUT" signal
EOF

echo
echo "=== Done. Outputs: ==="
echo "  $COMBINE_DIR/T2tt_mStop-*_mLSP-*/  (datacards, combine logs, shapes ROOT files)"
echo "  $LIMITS_CSV"
echo "  ${COMPARISON_OUT}.png / .pdf  (updated t2tt comparison plot, all mLSP points, like t2tt_equal_luminosity_comparison.png)"
echo "  ${COMPARISON_OUT}_values.csv"
echo "  $SEARCHBIN_COMBINED_OUT  (primary 174-bin plot: all mLSP signals summed into one 'signal' line)"
echo "  $SEARCHBIN_PLOT_DIR/T2tt_mStop-*_mLSP-*.png  (one 174-bin plot per mass point)"
