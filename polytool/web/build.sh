#!/bin/bash
set -e
cd "$(dirname "$0")"
~/.cargo/bin/wasm-pack build --target no-modules --out-dir pkg
echo ""
echo "Build complete ($(du -h pkg/*.wasm | cut -f1) WASM binary)."
echo "To test: cd web && python3 -m http.server"
