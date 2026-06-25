#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
package_dir="$(CDPATH= cd -- "$script_dir/../.." && pwd)"
manifest="$script_dir/manifest.tsv"

mkdir -p "$script_dir/json"

while IFS=$'\t' read -r slug _title _features args _recurrence rows_file json_file; do
  if [[ "$slug" == "slug" ]]; then
    continue
  fi

  # The manifest stores CLI arguments exactly as they should be passed.
  # shellcheck disable=SC2086
  cargo run -q --manifest-path "$package_dir/Cargo.toml" --bin polytool -- \
    recurrence --json $args \
    < "$script_dir/$rows_file" \
    > "$script_dir/$json_file"

  echo "wrote $json_file"
done < "$manifest"
