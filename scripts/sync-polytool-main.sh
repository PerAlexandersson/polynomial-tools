#!/usr/bin/env bash
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

if [[ $(git branch --show-current) != "master" ]]; then
    echo "error: run this script from the canonical master branch" >&2
    exit 2
fi

if [[ -n $(git status --porcelain -- polytool) ]]; then
    echo "error: polytool/ has uncommitted changes" >&2
    exit 2
fi

git fetch origin master main

if [[ $(git rev-parse master) != $(git rev-parse origin/master) ]]; then
    echo "error: local master and origin/master differ; push or update master first" >&2
    exit 2
fi

split_commit=$(git subtree split --prefix=polytool master)
git push origin "$split_commit:refs/heads/main"
git fetch origin main:refs/remotes/origin/main

echo "polytool main synchronized at $split_commit"
