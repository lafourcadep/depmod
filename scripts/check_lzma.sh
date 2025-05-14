#!/bin/bash

set -e

COMPILER=$1

if [ -z $COMPILER ]; then
    COMPILER="c++"
fi

CODE=$(cat <<EOF
#include <lzma.h>
int main() {
    lzma_stream strm = LZMA_STREAM_INIT;
    return lzma_easy_encoder(&strm, 6, LZMA_CHECK_CRC64);
}
EOF
)

TMPDIR=$(mktemp -d)
echo "$CODE" > "$TMPDIR/test.cpp"

if $COMPILER "$TMPDIR/test.cpp" -llzma -o "$TMPDIR/test" >/dev/null 2>&1; then
    rm -rf "$TMPDIR"
    exit 0
else
    rm -rf "$TMPDIR"
    exit 1
fi
