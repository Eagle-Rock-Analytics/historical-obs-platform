#!/bin/bash
# Run stnlist_update_merge.py for all networks.
# Usage: ./run_stnlist_update_merge.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

NETWORKS=(
    ASOSAWOS CAHYDRO CIMIS CW3E CDEC CNRFC CRN CWOP HADS
    HNXWFO HOLFUY HPWREN LOXWFO MAP MTRWFO NCAWOS NOS-NWLON
    NOS-PORTS otherisd RAWS SGXWFO SHASAVAL VCAPCD MARITIME
    NDBC SCAN SNOTEL VALLEYWATER
)

for network in "${NETWORKS[@]}"; do
    echo "Updating merge station list for $network..."
    python "$SCRIPT_DIR/stnlist_update_merge.py" "$network"
done

echo "Done."
