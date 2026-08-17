#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <repository-root> <output.tsv>" >&2
  exit 2
fi

repo_root="$(cd "$1" && pwd)"
output="$2"
mkdir -p "$(dirname "$output")"

git_commit="unavailable"
git_state="not_a_git_worktree"
if git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  if git -C "$repo_root" rev-parse --verify HEAD >/dev/null 2>&1; then
    git_commit="$(git -C "$repo_root" rev-parse HEAD)"
  else
    git_commit="unborn"
  fi
  if [[ -n "$(git -C "$repo_root" status --porcelain --untracked-files=all)" ]]; then
    git_state="modified_or_untracked"
  else
    git_state="clean"
  fi
fi

source_sha="$(
  cd "$repo_root"
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    source_list_command=(git ls-files --cached --others --exclude-standard -z)
  else
    source_list_command=(find . -type f ! -path '*/__pycache__/*' ! -path '*/output/*' -print0)
  fi
  while IFS= read -r -d '' source_file; do
    printf '%s\t%s\n' "$source_file" "$(sha256sum "$source_file" | awk '{print $1}')"
  done < <("${source_list_command[@]}" | LC_ALL=C sort -z) | sha256sum | awk '{print $1}'
)"

tmp="${output}.tmp.$$"
{
  printf 'field\tvalue\n'
  printf 'git_commit\t%s\n' "$git_commit"
  printf 'git_worktree_state\t%s\n' "$git_state"
  printf 'source_snapshot_sha256\t%s\n' "$source_sha"
} > "$tmp"
mv "$tmp" "$output"
