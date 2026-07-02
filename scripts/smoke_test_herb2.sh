#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-data/raw/test_downloads/herb2}"
mkdir -p "$OUT_DIR"

DOWNLOAD_BASE="http://47.92.70.12/download/file/"
BROWSE_API="http://herb.ac.cn/chedi/api/"
REMOTE_ROOT="/www/wwwroot/47.92.70.12/HERB_web/static/download_data/V2"
NO_PROXY_LIST="47.92.70.12,herb.ac.cn"
CURL=(curl -L -sS --fail --connect-timeout 20 --max-time 120 --retry 2 --noproxy "$NO_PROXY_LIST")

download_file() {
  local remote_name="$1"
  local output_name="$2"
  local output_path="$OUT_DIR/$output_name"
  "${CURL[@]}" --get "$DOWNLOAD_BASE" \
    --data-urlencode "file_path=$REMOTE_ROOT/$remote_name" \
    -o "$output_path"
  test -s "$output_path"
  awk -F'\t' 'NR==1 { if (NF < 5) exit 1 }' "$output_path"
  printf "%s\t%s\t%s\n" "$output_path" "$(wc -c < "$output_path" | tr -d ' ')" "$(md5 -q "$output_path")"
}

download_file "HERB_herb_info_v2.txt" "HERB_herb_info_v2.txt"
download_file "HERB_formula_info_v2.txt" "HERB_formula_info_v2.txt"

"${CURL[@]}" "$BROWSE_API" \
  -H 'Content-Type: application/json' \
  --data '{"label":"Herb","page":1,"page_size":15,"func_name":"browse_api"}' \
  -o "$OUT_DIR/herb_browse_page1.json"

test -s "$OUT_DIR/herb_browse_page1.json"
grep -q '"HERB' "$OUT_DIR/herb_browse_page1.json"
printf "%s\t%s\t%s\n" "$OUT_DIR/herb_browse_page1.json" "$(wc -c < "$OUT_DIR/herb_browse_page1.json" | tr -d ' ')" "$(md5 -q "$OUT_DIR/herb_browse_page1.json")"
