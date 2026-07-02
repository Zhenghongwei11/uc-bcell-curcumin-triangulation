#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-data/raw/test_downloads/etcm}"
mkdir -p "$OUT_DIR"

BROWSE_API="http://www.tcmip.cn:18124/home/browse/"
DETAIL_API="http://www.tcmip.cn:18124/home/detail/"
NO_PROXY_LIST="www.tcmip.cn"
CURL=(curl -L -sS --fail --connect-timeout 20 --max-time 120 --retry 2 --noproxy "$NO_PROXY_LIST")

browse_payload() {
  local output_name="$1"
  local type_name="$2"
  local page_size="${3:-20}"
  local output_path="$OUT_DIR/$output_name"
  "${CURL[@]}" "$BROWSE_API" \
    -H 'Content-Type: application/json' \
    --data "{\"type\":\"$type_name\",\"page\":1,\"pageSize\":$page_size,\"language\":\"en\"}" \
    -o "$output_path"
  test -s "$output_path"
  grep -q '"success' "$output_path"
  grep -q "\"type\":\"$type_name\"\\|\"type\": \"$type_name\"" "$output_path"
  printf "%s\t%s\t%s\n" "$output_path" "$(wc -c < "$output_path" | tr -d ' ')" "$(md5 -q "$output_path")"
}

browse_payload "formula_page1.json" "traditional_chinese_medicine_formula" 20
browse_payload "herb_page1.json" "herb" 20
browse_payload "ingredient_page1.json" "ingredient" 20

DETAIL_ID="${ETCM2_FORMULA_DETAIL_ID:-WuWeiZiWan10}"

"${CURL[@]}" --get "$DETAIL_API" \
  --data-urlencode "id=$DETAIL_ID" \
  --data-urlencode 'type=traditional_chinese_medicine_formula' \
  --data-urlencode 'language=en' \
  -o "$OUT_DIR/${DETAIL_ID}_detail.json"

test -s "$OUT_DIR/${DETAIL_ID}_detail.json"
grep -q '"code": 1' "$OUT_DIR/${DETAIL_ID}_detail.json"
grep -q '"base_information"' "$OUT_DIR/${DETAIL_ID}_detail.json"
grep -q '"ingredient_target"' "$OUT_DIR/${DETAIL_ID}_detail.json"
printf "%s\t%s\t%s\n" "$OUT_DIR/${DETAIL_ID}_detail.json" "$(wc -c < "$OUT_DIR/${DETAIL_ID}_detail.json" | tr -d ' ')" "$(md5 -q "$OUT_DIR/${DETAIL_ID}_detail.json")"
