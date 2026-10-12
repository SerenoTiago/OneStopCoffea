#!/usr/bin/env bash
# End-to-end Run-3 2024 T2tt luminosity projection at 137 fb^-1:
#   .result files -> datacards -> combine -> comparison plot + 174-bin plot
#
# Run this from the repository root in the normal analysis environment:
#   scripts/combine/run_full_pipeline.sh
#
# The CMS-SUS-19-006 Run-2 result is used only as an external benchmark. The
# local inputs, normalization, labels, and output names remain explicitly Run 3.
#
# Optional overrides (env vars):
#   SIGNAL_DIR      default: analysis_products/results/2026-10-11_run3_2024_t2tt_grid_expanded_signals
#                            (--test: analysis_products/old_results/26-09-10_hadsusy_signal_10k)
#   BACKGROUND_DIR  default: analysis_products/production_runs/2026-10-07_012944_hadronic_susy/background_results
#                            (complete 2024 backgrounds, same selection as the signals)
#                            (--test: analysis_products/old_results/26-09-10_hadsusy_backgrounds_10k)
#   BACKGROUND_GLOB default: *.result
#   SIGNAL_GLOB     default: *.result
#   COMBINE_DIR     default: analysis_products/combine/2026-10-11_run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb
#   TARGET_LUMINOSITY_FB default: 137.0
#   JOBS            default: 4 concurrent mass points
#   LOCAL_LUMINOSITY_FB  deprecated alias for TARGET_LUMINOSITY_FB
#   COLLISION_ENERGY_TEV default: 13.6
#   SEARCHBIN_PLOT_DIR   default: $COMBINE_DIR/search_bin_diagnostics
#   SEARCHBIN_COMBINED_OUT default: $COMBINE_DIR/search_bin_diagnostics_combined.png
#                          (primary plot: all mLSP signals summed into one "signal" line)
#   ANALYSIS_PYTHON default: $PWD/.venv-local10k/bin/python

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
    DEFAULT_SIGNAL_DIR="analysis_products/results/2026-10-11_run3_2024_t2tt_grid_expanded_signals"
    DEFAULT_BACKGROUND_DIR="analysis_products/production_runs/2026-10-07_012944_hadronic_susy/background_results"
fi

SIGNAL_DIR="${SIGNAL_DIR:-$DEFAULT_SIGNAL_DIR}"
BACKGROUND_DIR="${BACKGROUND_DIR:-$DEFAULT_BACKGROUND_DIR}"
SIGNAL_GLOB="${SIGNAL_GLOB:-*.result}"
BACKGROUND_GLOB="${BACKGROUND_GLOB:-*.result}"
MERGED_DIR="${MERGED_DIR:-.application_data/combine_inputs_2026-10-11_t2tt_grid_expanded_1007bkg}"
CONFIG_DIR="${CONFIG_DIR:-.application_data/generated_postprocess/2026-10-11_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb}"
COMBINE_DIR="${COMBINE_DIR:-analysis_products/combine/2026-10-11_run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb}"
LIMITS_CSV="${LIMITS_CSV:-$COMBINE_DIR/limits_expected.csv}"
SEARCHBIN_PLOT_DIR="${SEARCHBIN_PLOT_DIR:-$COMBINE_DIR/search_bin_diagnostics}"
SEARCHBIN_COMBINED_OUT="${SEARCHBIN_COMBINED_OUT:-$COMBINE_DIR/search_bin_diagnostics_combined.png}"
TARGET_LUMINOSITY_FB="${TARGET_LUMINOSITY_FB:-${LOCAL_LUMINOSITY_FB:-137.0}}"
COLLISION_ENERGY_TEV="${COLLISION_ENERGY_TEV:-13.6}"
JOBS="${JOBS:-4}"
ANALYSIS_PYTHON="${ANALYSIS_PYTHON:-$PWD/.venv-local10k/bin/python}"
ALLOW_MISSING_FLAG=""
if [[ "$TEST_MODE" == "true" ]]; then
    ALLOW_MISSING_FLAG="--allow-missing-backgrounds"
fi

echo "Mode: $([[ "$TEST_MODE" == "true" ]] && echo '--test (10k local)' || echo 'full condor rerun')"
echo "SIGNAL_DIR=$SIGNAL_DIR"
echo "BACKGROUND_DIR=$BACKGROUND_DIR"
echo "SIGNAL_GLOB=$SIGNAL_GLOB"
echo "BACKGROUND_GLOB=$BACKGROUND_GLOB"
echo "Run-3 luminosity projection: $TARGET_LUMINOSITY_FB fb^-1 at $COLLISION_ENERGY_TEV TeV"

echo "=== 0/4: merging signal + background result files into $MERGED_DIR ==="
if [[ "$MERGED_DIR" != .application_data/* ]]; then
    echo "ERROR: MERGED_DIR must remain under .application_data/: $MERGED_DIR" >&2
    exit 2
fi
mkdir -p "$MERGED_DIR"
if find "$MERGED_DIR" -mindepth 1 -maxdepth 1 ! -type l -print -quit | grep -q .; then
    echo "ERROR: refusing to clear non-symlink content from $MERGED_DIR" >&2
    exit 2
fi
find "$MERGED_DIR" -mindepth 1 -maxdepth 1 -type l -delete
ln -sf "$(realpath "$SIGNAL_DIR")"/$SIGNAL_GLOB "$MERGED_DIR"/
ln -sf "$(realpath "$BACKGROUND_DIR")"/$BACKGROUND_GLOB "$MERGED_DIR"/

echo "=== 1/4: generating + running per-signal datacards ==="
"$ANALYSIS_PYTHON" scripts/combine/make_split_signal_datacards.py \
    "$MERGED_DIR" --output-dir "$COMBINE_DIR" --config-dir "$CONFIG_DIR" \
    --run --overwrite --luminosity-override "$TARGET_LUMINOSITY_FB" \
    --jobs "$JOBS" $ALLOW_MISSING_FLAG

echo "=== 2/4: running combine (host apptainer) ==="
python3 scripts/combine/run_split_signal_combine.py \
    "$COMBINE_DIR" --keep-going --expected-only --jobs "$JOBS"

echo "=== 3/4: collecting limits + plotting the Run-3 expected-limit surface ==="
"$ANALYSIS_PYTHON" scripts/combine/collect_split_signal_limits.py \
    "$COMBINE_DIR" --output "$LIMITS_CSV" --expected-only
"$ANALYSIS_PYTHON" scripts/combine/plot_split_signal_limits.py \
    "$LIMITS_CSV" --output "$COMBINE_DIR/expected_limits_2d" \
    --quantity exp --paper-scale

echo "=== 4/4: 174-bin Run-3 search-region diagnostic plots ==="
mkdir -p "$SEARCHBIN_PLOT_DIR" "$(dirname "$SEARCHBIN_COMBINED_OUT")"

plot_diag() {
    local shapes_file="$1" output="$2" signal_hist="$3"
    local background_args
    background_args=$("$ANALYSIS_PYTHON" -c "
import uproot
keys = {k.split(';')[0] for k in uproot.open('$shapes_file').keys()}
defaults = ['qcd_inclusive_2024','tt_hadronic_2024','tt_semileptonic_2024','wjets_2024','zjets_2024']
print(' '.join(f'--background {b}' for b in defaults if b in keys))
")
    "$ANALYSIS_PYTHON" scripts/combine/plot_search_bin_diagnostics.py \
        "$shapes_file" --signal "$signal_hist" $background_args -o "$output" \
        --lumi "$TARGET_LUMINOSITY_FB" --energy "$COLLISION_ENERGY_TEV"
}

# One plot per mass point.
shapes_files=()
for point_dir in "$COMBINE_DIR"/mStop-*/T2tt_mStop-*_mLSP-*; do
    point=$(basename "$point_dir")
    shapes_file="$point_dir/shapes_bins_hadsusy.root"
    [[ -f "$shapes_file" ]] || continue
    shapes_files+=("$shapes_file")
    signal_hist=$("$ANALYSIS_PYTHON" -c "
import uproot
keys = {k.split(';')[0] for k in uproot.open('$shapes_file').keys()}
print(next(k for k in keys if k.startswith('2stop_')))
")
    plot_diag "$shapes_file" "$SEARCHBIN_PLOT_DIR/${point}.png" "$signal_hist"
done

# Primary combined plot: all mLSP signals summed into one "signal" line.
COMBINED_SHAPES="$COMBINE_DIR/shapes_all_signals_combined.root"
"$ANALYSIS_PYTHON" scripts/combine/build_combined_signal_shapes.py \
    "${shapes_files[@]}" -o "$COMBINED_SHAPES"
plot_diag "$COMBINED_SHAPES" "$SEARCHBIN_COMBINED_OUT" signal
echo
echo "=== Done. Outputs: ==="
echo "  $COMBINE_DIR/mStop-*/T2tt_mStop-*_mLSP-*/  (datacards, combine logs, shapes ROOT files)"
echo "  $LIMITS_CSV"
echo "  $COMBINE_DIR/expected_limits_2d-paper_scale.png / .pdf"
echo "  $COMBINE_DIR/expected_limits_2d-paper_scale_values.csv"
echo "  $SEARCHBIN_COMBINED_OUT  (primary 174-bin plot: all mLSP signals summed into one 'signal' line)"
echo "  $SEARCHBIN_PLOT_DIR/T2tt_mStop-*_mLSP-*.png  (one 174-bin plot per mass point)"
