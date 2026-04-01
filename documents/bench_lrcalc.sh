#!/bin/bash
# Benchmark: Ehrhart polynomial computation via lrcalc vs our kostka tool.
#
# For each test case (lambda, weight), compute the Ehrhart polynomial by:
# 1. Determining dimension d
# 2. Evaluating K(n*lambda, n*weight) at n=1,...,d+1 using lrcalc
#    (via the Kostka-to-LR reduction)
# 3. Timing the total lrcalc evaluations
#
# Then compare against `kostka ehrhart` timing.

set -e

KOSTKA_BIN="${KOSTKA_BIN:-cargo run --release -p kostka --}"

# Kostka-to-LR reduction: K_{lambda, mu} = c^rho_{lambda, sigma}
# where rho/sigma = disjoint union of single rows (mu_1), (mu_2), ...
#
# For weight mu = (w_1, w_2, ..., w_k):
#   widths = (w_1, w_2, ..., w_k)
#   col_offsets[i] = sum of widths[i+1..k]
#   rho row i = col_offsets[i] + w_i
#   sigma row i = col_offsets[i]  (or 0 if col_offset is 0 for last shape)
#
# For straight Kostka K_{lambda, mu} (no nu):
# shapes = single rows (w_1), (w_2), ..., (w_k)

kostka_via_lrcalc() {
    local lambda="$1"  # space-separated parts of n*lambda
    local weight="$2"  # space-separated parts of n*weight

    # Build rho and sigma from the weight parts
    local -a w_arr=($weight)
    local k=${#w_arr[@]}

    # Column offsets: col_offsets[i] = sum of w[i+1..k-1]
    local -a col_off
    col_off[$((k-1))]=0
    for ((i=k-2; i>=0; i--)); do
        col_off[$i]=$((col_off[$((i+1))] + w_arr[$((i+1))]))
    done

    # rho: row i has width col_off[i] + w_arr[i]
    # sigma: row i has width col_off[i] (0 for last row)
    local rho_parts=""
    local sigma_parts=""
    for ((i=0; i<k; i++)); do
        local rho_i=$((col_off[$i] + w_arr[$i]))
        local sig_i=${col_off[$i]}
        rho_parts="$rho_parts $rho_i"
        if [ $sig_i -gt 0 ]; then
            sigma_parts="$sigma_parts $sig_i"
        fi
    done
    rho_parts=$(echo $rho_parts | xargs)
    sigma_parts=$(echo $sigma_parts | xargs)

    # c^rho_{lambda, sigma}
    if [ -z "$sigma_parts" ]; then
        # sigma is empty: c^rho_{lambda, ∅} = delta_{rho, lambda}
        # (This shouldn't happen for interesting cases)
        echo "1"
    else
        lrcalc lrcoef $rho_parts - $lambda - $sigma_parts 2>/dev/null
    fi
}

scale_partition() {
    local n=$1
    shift
    local parts="$@"
    local result=""
    for p in $parts; do
        result="$result $((n * p))"
    done
    echo $result | xargs
}

echo "=== Benchmark: lrcalc vs kostka for Ehrhart polynomial ==="
echo ""

# Test cases: lambda weight
declare -a CASES=(
    "3 2 1|1 1 1 1 1 1"
    "3 3 3|1 1 1 1 1 1 1 1 1"
    "4 3 2|1 1 1 1 1 1 1 1 1"
    "4 2 2 2|2 1 1 1 1 1 1 1 1"
)

for case in "${CASES[@]}"; do
    IFS='|' read -r LAM_STR WT_STR <<< "$case"
    LAM_COMMA=$(echo $LAM_STR | tr ' ' ',')
    WT_COMMA=$(echo $WT_STR | tr ' ' ',')

    echo "--- lambda=($LAM_COMMA) weight=($WT_COMMA) ---"

    # Get dimension from kostka
    DIM=$($KOSTKA_BIN ehrhart --lambda "$LAM_COMMA" -w "$WT_COMMA" 2>/dev/null | grep "degree:" | awk '{print $2}')
    echo "  dimension: $DIM"
    NPOINTS=$((DIM + 1))

    # Time lrcalc evaluations
    echo -n "  lrcalc ($NPOINTS evals): "
    T0=$(date +%s%N)
    for ((n=1; n<=NPOINTS; n++)); do
        NLAM=$(scale_partition $n $LAM_STR)
        NWT=$(scale_partition $n $WT_STR)
        K=$(kostka_via_lrcalc "$NLAM" "$NWT")
    done
    T1=$(date +%s%N)
    LRCALC_MS=$(( (T1 - T0) / 1000000 ))
    echo "${LRCALC_MS} ms"

    # Time kostka ehrhart (with reciprocity)
    echo -n "  kostka ehrhart: "
    T0=$(date +%s%N)
    $KOSTKA_BIN ehrhart --lambda "$LAM_COMMA" -w "$WT_COMMA" > /dev/null 2>&1
    T1=$(date +%s%N)
    KOSTKA_MS=$(( (T1 - T0) / 1000000 ))
    echo "${KOSTKA_MS} ms"

    if [ $KOSTKA_MS -gt 0 ]; then
        echo "  speedup: ~$((LRCALC_MS / KOSTKA_MS))x"
    fi
    echo ""
done
