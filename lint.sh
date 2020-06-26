#!/bin/sh

SRC_DIR=source
INC_DIR=include
LINT_DIR=lint
LUA_INC_DIR="$LINT_DIR/lua"
NAUTY_INC_DIR="$LINT_DIR/nauty"

# clang-tidy
echo "running clang-tidy..."

run_clang_tidy() {
  echo "parsing $1..."

  clang-tidy -checks=-*,modernize-*,-modernize-use-trailing-return-type \
    "$SRC_DIR/$1" -header-filter="$INC_DIR/.*" \
    -- "$SRC_DIR/$1" \
    -I "$INC_DIR" -isystem "$LUA_INC_DIR" -isystem "$NAUTY_INC_DIR" \
    > "$LINT_DIR/clang-tidy-${1%.cpp}.txt" 2> /dev/null
}

for f in "$SRC_DIR"/*.cpp; do
  run_clang_tidy "$(basename "$f")"
done

# cppcheck
echo "running cppcheck..."

cppcheck --quiet \
  --std=c++11 \
  --enable=all --suppress=missingIncludeSystem \
  -ULUA_USER_H -UNAUTY_IN_MAGMA --force \
  -I "$INC_DIR"  -I "$LUA_INC_DIR" -I "$NAUTY_INC_DIR" \
  "$SRC_DIR" 2> "$LINT_DIR"/cppcheck.txt

# cppclean
echo "running cppclean..."

cppclean \
  -i "$INC_DIR" -s "$LUA_INC_DIR" -s "$NAUTY_INC_DIR" \
  "$INC_DIR" "$SRC_DIR" | grep -v "static data" > "$LINT_DIR/cppclean.txt"
