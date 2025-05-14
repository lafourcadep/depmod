#!/bin/bash

set -e

COMPILER=$1

if [ -z $COMPILER ]; then
    COMPILER="c++"
fi

CODE=$(cat <<EOF
#include <zlib.h>
int main() {
    return compress2(0, 0, 0, 0, 0);
}
EOF
)

TMPDIR=$(mktemp -d)
echo "$CODE" > "$TMPDIR/test.cpp"

if $COMPILER "$TMPDIR/test.cpp" -lz -o "$TMPDIR/test" >/dev/null 2>&1; then
    rm -rf "$TMPDIR"
    exit 0
else
    rm -rf "$TMPDIR"
    exit 1
fi
