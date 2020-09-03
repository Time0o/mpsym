#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage: $(basename $0) RESOURCE OUTFILE SYMBOL" 1>&2
fi

RESOURCE="$1"
OUTFILE="$2"
SYMBOL="$3"

DATA="$(hexdump -v -e '1/1 "0x%X,"' <<< "$(cat "$RESOURCE")")"

echo "#include <stddef.h>" >> "$OUTFILE"
echo "char const $SYMBOL[] = {$DATA};" >> "$OUTFILE"
echo "size_t const ${SYMBOL}_len = sizeof($SYMBOL);" >> "$OUTFILE"
