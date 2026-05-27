#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT_DIR/target/release/ktt-search"
LOG_DIR="${LOG_DIR:-$SCRIPT_DIR/logs}"

# Conservative defaults so SSH/Tailscale stays responsive.
JOBS="${JOBS:-1}"
MAX_JOBS="${MAX_JOBS:-2}"
SLICE_SECONDS="${SLICE_SECONDS:-60}"
SLEEP_SECONDS="${SLEEP_SECONDS:-20}"
STAGGER_SECONDS="${STAGGER_SECONDS:-30}"
NICE_N="${NICE_N:-15}"
MAX_STATES="${MAX_STATES:-200000}"
TAIL_LINES="${TAIL_LINES:-8}"
BUILD="${BUILD:-1}"
FLAG_MUTATIONS="${FLAG_MUTATIONS:-0}"

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || (( JOBS < 1 )); then
  echo "JOBS must be a positive integer" >&2
  exit 2
fi

if (( JOBS > MAX_JOBS )) && [[ "${ALLOW_MORE_JOBS:-0}" != "1" ]]; then
  echo "Refusing JOBS=$JOBS because MAX_JOBS=$MAX_JOBS." >&2
  echo "Use MAX_JOBS=$JOBS or ALLOW_MORE_JOBS=1 if you really want that many." >&2
  exit 2
fi

mkdir -p "$LOG_DIR"

if [[ "$BUILD" == "1" || ! -x "$BIN" ]]; then
  echo "[setup] building release binary"
  (cd "$ROOT_DIR" && timeout 60s nice -n 10 cargo build --release -p ktt-search)
fi

base_args=(--max-states "$MAX_STATES")
if [[ "$FLAG_MUTATIONS" != "1" ]]; then
  base_args+=(--no-flag-mutations)
fi
base_args+=("$@")

ionice_cmd=()
if command -v ionice >/dev/null 2>&1; then
  ionice_cmd=(ionice -c 3)
fi

pids=()

stop_children() {
  trap - INT TERM
  if ((${#pids[@]})); then
    echo
    echo "[stop] stopping workers: ${pids[*]}" >&2
    kill "${pids[@]}" 2>/dev/null || true
    pids=()
  fi
}
trap stop_children INT TERM

run_worker() {
  local worker_id="$1"
  local initial_delay=$(( (worker_id - 1) * STAGGER_SECONDS ))
  if (( initial_delay > 0 )); then
    echo "[worker $worker_id] stagger sleep ${initial_delay}s"
    sleep "$initial_delay"
  fi

  while true; do
    local stamp log code
    stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    log="$LOG_DIR/ktt-search-w${worker_id}-${stamp}.log"

    echo "[worker $worker_id] start slice ${SLICE_SECONDS}s -> $log"
    set +e
    timeout "${SLICE_SECONDS}s" nice -n "$NICE_N" "${ionice_cmd[@]}" "$BIN" "${base_args[@]}" >"$log" 2>&1
    code=$?
    set -e

    if grep -q "NEGATIVE COEFFICIENT FOUND" "$log"; then
      echo "[worker $worker_id] NEGATIVE COEFFICIENT FOUND; see $log"
      exit 42
    fi

    case "$code" in
      0)
        echo "[worker $worker_id] search completed; see $log"
        tail -n "$TAIL_LINES" "$log" || true
        exit 0
        ;;
      124)
        echo "[worker $worker_id] slice timeout; last $TAIL_LINES lines:"
        tail -n "$TAIL_LINES" "$log" || true
        ;;
      *)
        echo "[worker $worker_id] command failed with exit $code; see $log" >&2
        tail -n "$TAIL_LINES" "$log" || true
        exit "$code"
        ;;
    esac

    echo "[worker $worker_id] sleep ${SLEEP_SECONDS}s"
    sleep "$SLEEP_SECONDS"
  done
}

echo "[config] JOBS=$JOBS SLICE_SECONDS=$SLICE_SECONDS SLEEP_SECONDS=$SLEEP_SECONDS NICE_N=$NICE_N MAX_STATES=$MAX_STATES FLAG_MUTATIONS=$FLAG_MUTATIONS"
echo "[config] extra args: ${*:-<none>}"

for worker_id in $(seq 1 "$JOBS"); do
  run_worker "$worker_id" &
  pids+=("$!")
done

if wait -n; then
  status=0
else
  status=$?
fi
stop_children
exit "$status"
