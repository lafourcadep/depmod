#!/bin/bash

set -e

COMPILER=$1

if [ -z $COMPILER ]; then
    COMPILER="c++"
fi

CODE=$(cat <<EOF
#include <bzlib.h>
int main() {
    bz_stream strm = {0};
    return BZ2_bzCompressInit(&strm, 9, 0, 0);
}
EOF
)

TMPDIR=$(mktemp -d)
echo "$CODE" > "$TMPDIR/test.cpp"

if $COMPILER "$TMPDIR/test.cpp" -lbz2 -o "$TMPDIR/test" >/dev/null 2>&1; then
    rm -rf "$TMPDIR"
    exit 0
else
    rm -rf "$TMPDIR"
    exit 1
fi
