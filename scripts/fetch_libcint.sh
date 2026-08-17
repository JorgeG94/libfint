#!/usr/bin/env bash
# Fetch libcint, for verification only.
#
#     scripts/fetch_libcint.sh [ref]      # default: master
#
# libfint contains no C.  The bit-identity suite needs libcint to compare
# against, so CI fetches it here rather than vendoring or submoduling it.
#
# The default is upstream's master rather than a pinned tag, deliberately.
# This suite is a drift detector: when it goes red, either upstream changed
# something the Fortran generator has not been told about, or upstream fixed
# something.  Either way that is news, and a pin would suppress it.  Pass a
# ref to reproduce an older run.
set -euo pipefail

REF="${1:-master}"
REPO="${LIBCINT_REPO:-https://github.com/sunqm/libcint}"
DEST="${LIBCINT_DIR:-extern/libcint}"

if [ -d "$DEST" ]; then
    echo "fetch_libcint: $DEST exists, leaving it alone"
    exit 0
fi

mkdir -p "$(dirname "$DEST")"

# A tarball, not a clone: no git history is wanted and the download is a
# tenth of the size.  Falls back to git for refs the archive endpoint will
# not serve.
if command -v curl >/dev/null 2>&1; then
    echo "fetch_libcint: $REPO @ $REF -> $DEST"
    tmp="$(mktemp -d)"
    if curl -fsSL "$REPO/archive/$REF.tar.gz" -o "$tmp/libcint.tar.gz"; then
        mkdir -p "$DEST"
        tar -xzf "$tmp/libcint.tar.gz" -C "$DEST" --strip-components=1
        rm -rf "$tmp"
        echo "fetch_libcint: ok"
        exit 0
    fi
    rm -rf "$tmp"
    echo "fetch_libcint: archive fetch failed, falling back to git"
fi

git clone --depth 1 --branch "$REF" "$REPO" "$DEST"
echo "fetch_libcint: ok"
