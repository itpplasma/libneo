#!/bin/bash
#> Script to dynamically fetch golden record odeint implementation
#> from the latest tagged version for comparison testing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GOLDEN_FILE="$SCRIPT_DIR/golden_odeint_dynamic.f90"

# Detect OS and set appropriate sed command
if [[ "$OSTYPE" == "darwin"* ]]; then
    SED_CMD="gsed"
else
    SED_CMD="sed"
fi

echo "Fetching golden record from latest tag..."

# Get the latest tag from remote (works even on shallow clones)
LATEST_TAG=$(git ls-remote --tags --sort=-v:refname origin | head -1 | sed 's|.*refs/tags/||')
echo "Latest tag: $LATEST_TAG"

# Fetch only that tagged commit (no-op if already available locally)
git fetch origin tag "$LATEST_TAG" --no-tags

# Fetch the old implementation
echo "Fetching odeint_allroutines.f90 from $LATEST_TAG..."
git show "$LATEST_TAG:src/odeint_allroutines.f90" > "$GOLDEN_FILE"

# Rename modules to avoid conflicts
echo "Adapting golden record for testing..."
$SED_CMD -i 's/module odeint_mod/module odeint_golden_mod/g' "$GOLDEN_FILE"
$SED_CMD -i 's/module odeint_allroutines_sub/module odeint_golden_sub/g' "$GOLDEN_FILE"
$SED_CMD -i 's/use odeint_mod/use odeint_golden_mod/g' "$GOLDEN_FILE"

echo "Golden record ready: $GOLDEN_FILE"
echo "Use 'make test_golden_record' to run comparison tests"