#!/usr/bin/env bash
# Run all T2tt jet-bin-fix diagnostic scripts in one shot.
#
# Enter the container first (venv/uv only exist there), then run from the repo root:
#   ./setup.sh
#   scripts/diagnostics/run_all.sh            # uses the full-stats condor rerun output
#   scripts/diagnostics/run_all.sh --test      # uses the quick local 10k-event result files
#
# Optional overrides (env vars, take precedence over --test defaults):
#   BEFORE_SIGNAL   default: analysis_products/old_results/26-08-30_hadsusy_1250_run2
#   AFTER_SIGNAL    default: analysis_products/old_results/26-09-10_hadsusy_jetbins_fixed/signal
#                   (--test: analysis_products/old_results/26-09-10_hadsusy_signal_10k)
#   BEFORE_BKG      default: analysis_products/old_results/26-08-30_hadsusy_1250_run2
#   AFTER_BKG       default: analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed
#                   (--test: analysis_products/old_results/26-09-10_hadsusy_backgrounds_10k)
#   OUT_DIR         default: analysis_products/plots/diagnostics

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

TEST_MODE=false
for arg in "$@"; do
    if [[ "$arg" == "--test" ]]; then
        TEST_MODE=true
    fi
done

export UV_PROJECT_ENVIRONMENT="${UV_PROJECT_ENVIRONMENT:-$PWD/.venv-local10k}"
export PYTHONPATH="$PWD"

if [[ "$TEST_MODE" == "true" ]]; then
    DEFAULT_AFTER_SIGNAL="analysis_products/old_results/26-09-10_hadsusy_signal_10k"
    DEFAULT_AFTER_BKG="analysis_products/old_results/26-09-10_hadsusy_backgrounds_10k"
else
    DEFAULT_AFTER_SIGNAL="analysis_products/old_results/26-09-10_hadsusy_jetbins_fixed/signal"
    DEFAULT_AFTER_BKG="analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed"
fi

BEFORE_SIGNAL="${BEFORE_SIGNAL:-analysis_products/old_results/26-08-30_hadsusy_1250_run2}"
AFTER_SIGNAL="${AFTER_SIGNAL:-$DEFAULT_AFTER_SIGNAL}"
BEFORE_BKG="${BEFORE_BKG:-analysis_products/old_results/26-08-30_hadsusy_1250_run2}"
AFTER_BKG="${AFTER_BKG:-$DEFAULT_AFTER_BKG}"
OUT_DIR="${OUT_DIR:-analysis_products/plots/diagnostics}"

echo "Mode: $([[ "$TEST_MODE" == "true" ]] && echo '--test (10k local)' || echo 'full condor rerun')"
echo "AFTER_SIGNAL=$AFTER_SIGNAL"
echo "AFTER_BKG=$AFTER_BKG"

n_after_bkg=$(ls "$AFTER_BKG"/*.result 2>/dev/null | wc -l)
if [[ "$n_after_bkg" -lt 27 ]]; then
    echo "WARNING: $AFTER_BKG has only $n_after_bkg/27 result files. Results will be incomplete/preliminary." >&2
fi

echo "=== 1/3: Signal JetN before/after (1250/100) ==="
uv run --no-sync python3 scripts/diagnostics/plot_jetn_before_after.py \
    --before "$BEFORE_SIGNAL" \
    --after "$AFTER_SIGNAL" \
    --sample-glob signal_T2tt_1250_100 \
    --hist-name JetN \
    --output "$OUT_DIR/jetn_signal_1250_100_before_after.png"

echo "=== 2/3: Background JetN before/after (summed over all 5 background datasets) ==="
uv run --no-sync python3 scripts/diagnostics/plot_jetn_before_after.py \
    --before "$BEFORE_BKG" \
    --after "$AFTER_BKG" \
    --sample-glob "" \
    --hist-name JetN \
    --output "$OUT_DIR/jetn_backgrounds_before_after.png"

echo "=== 3/3: SR background projections vs CMS pre-fit background ==="
uv run --no-sync python3 scripts/diagnostics/compare_background_projections.py \
    --results "$AFTER_BKG" \
    --output "$OUT_DIR/local_vs_cms_projections"

echo
echo "=== Done. Outputs: ==="
echo "  $OUT_DIR/jetn_signal_1250_100_before_after.png"
echo "  $OUT_DIR/jetn_backgrounds_before_after.png"
echo "  $OUT_DIR/local_vs_cms_projections/{njet,nbjet,htmiss}_components_{absolute,shape}.png"
echo "  $OUT_DIR/local_vs_cms_projections/component_projection_values.csv"
echo "  $OUT_DIR/local_vs_cms_projections/{njet,nbjet,htmiss}_{absolute,shape}.png"
echo "  $OUT_DIR/local_vs_cms_projections/projection_values.csv"
