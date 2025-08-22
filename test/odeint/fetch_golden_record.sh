#!/bin/bash
#> Script to dynamically fetch golden record odeint implementation
#> from the latest tagged version for comparison testing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GOLDEN_FILE="$SCRIPT_DIR/golden_odeint_dynamic.f90"

echo "🔍 Fetching golden record from latest tag..."

# Get the latest tag
LATEST_TAG=$(git tag --sort=-v:refname | head -1)
echo "📌 Latest tag: $LATEST_TAG"

# Fetch the old implementation
echo "📥 Fetching odeint_allroutines.f90 from $LATEST_TAG..."
git show "$LATEST_TAG:src/odeint_allroutines.f90" > "$GOLDEN_FILE"

# Rename modules to avoid conflicts
echo "🔧 Adapting golden record for testing..."
sed -i 's/module odeint_mod/module odeint_golden_mod/g' "$GOLDEN_FILE"
sed -i 's/module odeint_allroutines_sub/module odeint_golden_sub/g' "$GOLDEN_FILE"
sed -i 's/use odeint_mod/use odeint_golden_mod/g' "$GOLDEN_FILE"

echo "✅ Golden record ready: $GOLDEN_FILE"
echo "🧪 Use 'make test_golden_record' to run comparison tests"